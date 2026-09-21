//! Read selection by an expression over SAM aux tags: the `--tag-filter`
//! grammar, its parser and an evaluator over any record type that can look up
//! a tag.
//!
//! The grammar is the aux-tag subset of the `samtools view -e` language:
//!
//! ```text
//! expr    := or
//! or      := and ( "||" and )*
//! and     := not ( "&&" not )*
//! not     := "!" not | cmp
//! cmp     := primary ( ( "==" | "!=" | "<" | "<=" | ">" | ">=" ) primary )?
//! primary := "(" expr ")" | "[" TAG "]" | "exists" "(" "[" TAG "]" ")" | NUMBER | STRING
//! ```
//!
//! A tag in a boolean position is true when the tag is present. Comparisons
//! are numeric between numbers and bytewise between strings; a comparison
//! with a missing tag is false, whatever the operator, as in samtools. A
//! comparison between a number and a string, or one involving an array, is
//! an error that names the tag, since a typo must not silently select or
//! drop every read.

use std::borrow::Cow;
use std::fmt;

/// A tag value as the evaluator sees it, independent of the record type.
#[derive(Debug, Clone, PartialEq)]
pub enum TagValue<'a> {
    /// The tag is not on the record.
    Missing,
    /// An integer of any width or a float (`i`, `f`).
    Num(f64),
    /// A string, hex string or single character (`Z`, `H`, `A`).
    Str(Cow<'a, [u8]>),
    /// A `B` array; only its presence can be tested.
    Array,
}

/// A record whose aux tags the filter can read.
pub trait TagSource {
    /// Returns the value of `tag`, or `Missing`. An unreadable tag is an error.
    fn tag_value(&self, tag: [u8; 2]) -> anyhow::Result<TagValue<'_>>;
}

/// A comparison operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Op {
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
}

impl Op {
    fn holds(self, ord: std::cmp::Ordering) -> bool {
        use std::cmp::Ordering::*;
        match self {
            Op::Eq => ord == Equal,
            Op::Ne => ord != Equal,
            Op::Lt => ord == Less,
            Op::Le => ord != Greater,
            Op::Gt => ord == Greater,
            Op::Ge => ord != Less,
        }
    }

    fn symbol(self) -> &'static str {
        match self {
            Op::Eq => "==",
            Op::Ne => "!=",
            Op::Lt => "<",
            Op::Le => "<=",
            Op::Gt => ">",
            Op::Ge => ">=",
        }
    }
}

/// A parsed expression.
#[derive(Debug, Clone, PartialEq)]
enum Expr {
    Tag([u8; 2]),
    Num(f64),
    Str(Vec<u8>),
    Not(Box<Expr>),
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    Cmp(Op, Box<Expr>, Box<Expr>),
}

/// One `--tag-filter` expression, compiled once and shared by every reader.
#[derive(Debug, Clone, PartialEq)]
pub struct TagFilter {
    text: String,
    expr: Expr,
}

/// The `--tag-filter` expressions of a run; a read is kept when every one
/// holds.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TagFilters {
    filters: Vec<TagFilter>,
}

impl TagFilters {
    /// Parses every expression; the first failure names its expression.
    pub fn parse(values: &[String]) -> anyhow::Result<Self> {
        let filters = values
            .iter()
            .map(|v| TagFilter::parse(v))
            .collect::<anyhow::Result<Vec<_>>>()?;
        Ok(TagFilters { filters })
    }

    /// Whether no expression was given.
    pub fn is_empty(&self) -> bool {
        self.filters.is_empty()
    }

    /// The expressions as written.
    pub fn texts(&self) -> impl Iterator<Item = &str> {
        self.filters.iter().map(|f| f.text.as_str())
    }

    /// Whether `record` satisfies every expression.
    pub fn keeps(&self, record: &impl TagSource) -> anyhow::Result<bool> {
        for f in &self.filters {
            if !f.keeps(record)? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

impl TagFilter {
    /// Parses one expression.
    pub fn parse(text: &str) -> anyhow::Result<Self> {
        let tokens = lex(text).map_err(|e| anyhow::anyhow!("--tag-filter {text:?}: {e}"))?;
        let mut parser = Parser { tokens, pos: 0 };
        let expr = parser
            .expr()
            .and_then(|e| {
                parser.expect_end()?;
                require_condition(&e)?;
                Ok(e)
            })
            .map_err(|e| anyhow::anyhow!("--tag-filter {text:?}: {e}"))?;
        Ok(TagFilter {
            text: text.to_string(),
            expr,
        })
    }

    /// The expression as written.
    pub fn text(&self) -> &str {
        &self.text
    }

    /// Whether `record` satisfies the expression.
    pub fn keeps(&self, record: &impl TagSource) -> anyhow::Result<bool> {
        eval_bool(&self.expr, record)
    }
}

/// A bare literal is a value, not a condition.
fn require_condition(e: &Expr) -> Result<(), String> {
    match e {
        Expr::Num(_) | Expr::Str(_) => Err("a literal is not a condition".to_string()),
        _ => Ok(()),
    }
}

/// A lexical token.
#[derive(Debug, Clone, PartialEq)]
enum Token {
    LParen,
    RParen,
    Tag([u8; 2]),
    Exists,
    Not,
    And,
    Or,
    Cmp(Op),
    Num(f64),
    Str(Vec<u8>),
}

impl fmt::Display for Token {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Token::LParen => write!(f, "("),
            Token::RParen => write!(f, ")"),
            Token::Tag(t) => write!(f, "[{}]", String::from_utf8_lossy(t)),
            Token::Exists => write!(f, "exists"),
            Token::Not => write!(f, "!"),
            Token::And => write!(f, "&&"),
            Token::Or => write!(f, "||"),
            Token::Cmp(op) => write!(f, "{}", op.symbol()),
            Token::Num(n) => write!(f, "{n}"),
            Token::Str(s) => write!(f, "\"{}\"", String::from_utf8_lossy(s)),
        }
    }
}

/// Splits `text` into tokens.
fn lex(text: &str) -> Result<Vec<Token>, String> {
    let b = text.as_bytes();
    let mut i = 0;
    let mut out = Vec::new();
    while i < b.len() {
        let c = b[i];
        match c {
            b' ' | b'\t' | b'\n' | b'\r' => i += 1,
            b'(' => {
                out.push(Token::LParen);
                i += 1;
            },
            b')' => {
                out.push(Token::RParen);
                i += 1;
            },
            b'[' => {
                let end = b[i..]
                    .iter()
                    .position(|&x| x == b']')
                    .map(|p| i + p)
                    .ok_or_else(|| format!("unterminated tag at position {}", i + 1))?;
                let name = &b[i + 1..end];
                if name.len() != 2
                    || !name[0].is_ascii_alphabetic()
                    || !name[1].is_ascii_alphanumeric()
                {
                    return Err(format!(
                        "invalid tag [{}] at position {} (a SAM tag is a letter followed by a letter or digit)",
                        String::from_utf8_lossy(name),
                        i + 1
                    ));
                }
                out.push(Token::Tag([name[0], name[1]]));
                i = end + 1;
            },
            b'!' if b.get(i + 1) == Some(&b'=') => {
                out.push(Token::Cmp(Op::Ne));
                i += 2;
            },
            b'!' => {
                out.push(Token::Not);
                i += 1;
            },
            b'=' if b.get(i + 1) == Some(&b'=') => {
                out.push(Token::Cmp(Op::Eq));
                i += 2;
            },
            b'<' if b.get(i + 1) == Some(&b'=') => {
                out.push(Token::Cmp(Op::Le));
                i += 2;
            },
            b'<' => {
                out.push(Token::Cmp(Op::Lt));
                i += 1;
            },
            b'>' if b.get(i + 1) == Some(&b'=') => {
                out.push(Token::Cmp(Op::Ge));
                i += 2;
            },
            b'>' => {
                out.push(Token::Cmp(Op::Gt));
                i += 1;
            },
            b'&' if b.get(i + 1) == Some(&b'&') => {
                out.push(Token::And);
                i += 2;
            },
            b'|' if b.get(i + 1) == Some(&b'|') => {
                out.push(Token::Or);
                i += 2;
            },
            b'"' => {
                let mut s = Vec::new();
                let mut j = i + 1;
                loop {
                    match b.get(j) {
                        None => return Err(format!("unterminated string at position {}", i + 1)),
                        Some(b'"') => break,
                        Some(b'\\') => match b.get(j + 1) {
                            Some(&e @ (b'"' | b'\\')) => {
                                s.push(e);
                                j += 2;
                            },
                            _ => {
                                return Err(format!(
                                    "invalid escape at position {} (only \\\" and \\\\ are accepted)",
                                    j + 1
                                ));
                            },
                        },
                        Some(&x) => {
                            s.push(x);
                            j += 1;
                        },
                    }
                }
                out.push(Token::Str(s));
                i = j + 1;
            },
            b'0'..=b'9' | b'-' | b'+' | b'.' => {
                let start = i;
                i += 1;
                while i < b.len() && matches!(b[i], b'0'..=b'9' | b'.' | b'e' | b'E' | b'-' | b'+')
                {
                    // A sign continues a number only directly after an exponent marker.
                    if matches!(b[i], b'-' | b'+') && !matches!(b[i - 1], b'e' | b'E') {
                        break;
                    }
                    i += 1;
                }
                let s = &text[start..i];
                let n: f64 = s
                    .parse()
                    .map_err(|_| format!("invalid number {s:?} at position {}", start + 1))?;
                out.push(Token::Num(n));
            },
            _ if c.is_ascii_alphabetic() => {
                let start = i;
                while i < b.len() && (b[i].is_ascii_alphanumeric() || b[i] == b'_') {
                    i += 1;
                }
                let word = &text[start..i];
                match word {
                    "exists" => out.push(Token::Exists),
                    _ => {
                        return Err(format!(
                            "unsupported word {word:?} at position {} (tags are written as [{word}]; \
                             functions other than exists and record fields are not supported)",
                            start + 1
                        ));
                    },
                }
            },
            _ => {
                return Err(format!(
                    "unexpected character {:?} at position {}",
                    c as char,
                    i + 1
                ));
            },
        }
    }
    Ok(out)
}

/// A recursive-descent parser over the token list.
struct Parser {
    tokens: Vec<Token>,
    pos: usize,
}

impl Parser {
    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.pos)
    }

    fn next(&mut self) -> Option<Token> {
        let t = self.tokens.get(self.pos).cloned();
        self.pos += 1;
        t
    }

    fn expect(&mut self, want: Token) -> Result<(), String> {
        match self.next() {
            Some(t) if t == want => Ok(()),
            Some(t) => Err(format!("expected {want} but found {t}")),
            None => Err(format!("expected {want} at the end of the expression")),
        }
    }

    fn expect_end(&mut self) -> Result<(), String> {
        match self.next() {
            None => Ok(()),
            Some(t) => Err(format!("unexpected {t} after the expression")),
        }
    }

    fn expr(&mut self) -> Result<Expr, String> {
        let mut left = self.and()?;
        while self.peek() == Some(&Token::Or) {
            self.pos += 1;
            let right = self.and()?;
            left = Expr::Or(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn and(&mut self) -> Result<Expr, String> {
        let mut left = self.not()?;
        while self.peek() == Some(&Token::And) {
            self.pos += 1;
            let right = self.not()?;
            left = Expr::And(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn not(&mut self) -> Result<Expr, String> {
        if self.peek() == Some(&Token::Not) {
            self.pos += 1;
            let inner = self.not()?;
            require_condition(&inner)?;
            return Ok(Expr::Not(Box::new(inner)));
        }
        self.cmp()
    }

    fn cmp(&mut self) -> Result<Expr, String> {
        let left = self.primary()?;
        if let Some(Token::Cmp(op)) = self.peek().cloned() {
            self.pos += 1;
            let right = self.primary()?;
            if matches!(left, Expr::Num(_) | Expr::Str(_))
                && matches!(right, Expr::Num(_) | Expr::Str(_))
            {
                return Err("a comparison needs a tag on at least one side".to_string());
            }
            return Ok(Expr::Cmp(op, Box::new(left), Box::new(right)));
        }
        Ok(left)
    }

    fn primary(&mut self) -> Result<Expr, String> {
        match self.next() {
            Some(Token::LParen) => {
                let e = self.expr()?;
                self.expect(Token::RParen)?;
                Ok(e)
            },
            Some(Token::Tag(t)) => Ok(Expr::Tag(t)),
            Some(Token::Exists) => {
                self.expect(Token::LParen)?;
                let tag = match self.next() {
                    Some(Token::Tag(t)) => t,
                    Some(t) => return Err(format!("exists() takes a tag, found {t}")),
                    None => return Err("exists() takes a tag".to_string()),
                };
                self.expect(Token::RParen)?;
                Ok(Expr::Tag(tag))
            },
            Some(Token::Num(n)) => Ok(Expr::Num(n)),
            Some(Token::Str(s)) => Ok(Expr::Str(s)),
            Some(t) => Err(format!("unexpected {t}")),
            None => Err("the expression ends early".to_string()),
        }
    }
}

/// Evaluates `e` as a condition.
fn eval_bool(e: &Expr, rec: &impl TagSource) -> anyhow::Result<bool> {
    Ok(match e {
        Expr::Tag(t) => rec.tag_value(*t)? != TagValue::Missing,
        Expr::Not(inner) => !eval_bool(inner, rec)?,
        Expr::And(a, b) => eval_bool(a, rec)? && eval_bool(b, rec)?,
        Expr::Or(a, b) => eval_bool(a, rec)? || eval_bool(b, rec)?,
        Expr::Cmp(op, a, b) => compare(*op, a, b, rec)?,
        Expr::Num(_) | Expr::Str(_) => {
            anyhow::bail!("a literal is not a condition")
        },
    })
}

/// Evaluates `e` as a value.
fn eval_value<'r>(
    e: &'r Expr,
    rec: &'r impl TagSource,
) -> anyhow::Result<(TagValue<'r>, Option<[u8; 2]>)> {
    Ok(match e {
        Expr::Tag(t) => (rec.tag_value(*t)?, Some(*t)),
        Expr::Num(n) => (TagValue::Num(*n), None),
        Expr::Str(s) => (TagValue::Str(Cow::Borrowed(s)), None),
        _ => anyhow::bail!("a condition cannot be compared with {}", op_context(e)),
    })
}

fn op_context(e: &Expr) -> &'static str {
    match e {
        Expr::Cmp(..) => "another comparison",
        Expr::Not(_) => "a negation",
        Expr::And(..) | Expr::Or(..) => "a boolean expression",
        Expr::Tag(_) | Expr::Num(_) | Expr::Str(_) => "a value",
    }
}

/// Compares two values. A missing tag makes the comparison false; a number
/// against a string, or an array, is an error naming the tag.
fn compare(op: Op, a: &Expr, b: &Expr, rec: &impl TagSource) -> anyhow::Result<bool> {
    let (left, ltag) = eval_value(a, rec)?;
    let (right, rtag) = eval_value(b, rec)?;
    let ord = match (&left, &right) {
        (TagValue::Missing, _) | (_, TagValue::Missing) => return Ok(false),
        (TagValue::Num(x), TagValue::Num(y)) => match x.partial_cmp(y) {
            Some(ord) => ord,
            None => return Ok(false),
        },
        (TagValue::Str(x), TagValue::Str(y)) => x.as_ref().cmp(y.as_ref()),
        _ => {
            let tag = ltag.or(rtag).unwrap_or(*b"??");
            anyhow::bail!(
                "--tag-filter compares [{}] ({}) {} {} ({}); compare a number with a number and a string with a string",
                String::from_utf8_lossy(&tag),
                kind(&left),
                op.symbol(),
                describe(&right),
                kind(&right)
            );
        },
    };
    Ok(op.holds(ord))
}

fn kind(v: &TagValue<'_>) -> &'static str {
    match v {
        TagValue::Missing => "missing",
        TagValue::Num(_) => "a number",
        TagValue::Str(_) => "a string",
        TagValue::Array => "an array",
    }
}

fn describe(v: &TagValue<'_>) -> String {
    match v {
        TagValue::Missing => "missing".to_string(),
        TagValue::Num(n) => n.to_string(),
        TagValue::Str(s) => format!("\"{}\"", String::from_utf8_lossy(s)),
        TagValue::Array => "an array".to_string(),
    }
}

/// Converts a borrowed BAM aux value. Strings are copied, since the borrowed
/// value does not outlive the lookup.
pub fn value_of_raw<'a>(
    value: &noodles_sam::alignment::record::data::field::Value<'_>,
) -> TagValue<'a> {
    use noodles_sam::alignment::record::data::field::Value as V;
    match value {
        V::Character(c) => TagValue::Str(Cow::Owned(vec![*c])),
        V::Int8(n) => TagValue::Num(f64::from(*n)),
        V::UInt8(n) => TagValue::Num(f64::from(*n)),
        V::Int16(n) => TagValue::Num(f64::from(*n)),
        V::UInt16(n) => TagValue::Num(f64::from(*n)),
        V::Int32(n) => TagValue::Num(f64::from(*n)),
        V::UInt32(n) => TagValue::Num(f64::from(*n)),
        V::Float(x) => TagValue::Num(f64::from(*x)),
        V::String(s) | V::Hex(s) => TagValue::Str(Cow::Owned(s.to_vec())),
        V::Array(_) => TagValue::Array,
    }
}

/// Converts an owned aux value.
pub fn value_of_buf(
    value: &noodles_sam::alignment::record_buf::data::field::Value,
) -> TagValue<'_> {
    use noodles_sam::alignment::record_buf::data::field::Value as V;
    match value {
        V::Character(c) => TagValue::Str(Cow::Owned(vec![*c])),
        V::Int8(n) => TagValue::Num(f64::from(*n)),
        V::UInt8(n) => TagValue::Num(f64::from(*n)),
        V::Int16(n) => TagValue::Num(f64::from(*n)),
        V::UInt16(n) => TagValue::Num(f64::from(*n)),
        V::Int32(n) => TagValue::Num(f64::from(*n)),
        V::UInt32(n) => TagValue::Num(f64::from(*n)),
        V::Float(x) => TagValue::Num(f64::from(*x)),
        V::String(s) | V::Hex(s) => TagValue::Str(Cow::Borrowed(s.as_ref())),
        V::Array(_) => TagValue::Array,
    }
}

/// A record without aux tags: every lookup is `Missing`.
pub struct NoTags;

impl TagSource for NoTags {
    fn tag_value(&self, _tag: [u8; 2]) -> anyhow::Result<TagValue<'_>> {
        Ok(TagValue::Missing)
    }
}

impl TagSource for noodles_bam::Record {
    fn tag_value(&self, tag: [u8; 2]) -> anyhow::Result<TagValue<'_>> {
        use noodles_sam::alignment::record::data::field::Tag;
        match self.data().get(&Tag::from(tag)) {
            None => Ok(TagValue::Missing),
            Some(Ok(value)) => Ok(value_of_raw(&value)),
            Some(Err(e)) => Err(anyhow::anyhow!(
                "tag {} could not be read: {e}",
                String::from_utf8_lossy(&tag)
            )),
        }
    }
}

impl TagSource for noodles_sam::alignment::record_buf::Data {
    fn tag_value(&self, tag: [u8; 2]) -> anyhow::Result<TagValue<'_>> {
        use noodles_sam::alignment::record::data::field::Tag;
        Ok(self
            .get(&Tag::from(tag))
            .map_or(TagValue::Missing, value_of_buf))
    }
}

impl TagSource for noodles_sam::alignment::RecordBuf {
    fn tag_value(&self, tag: [u8; 2]) -> anyhow::Result<TagValue<'_>> {
        self.data().tag_value(tag)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    struct Rec(HashMap<[u8; 2], TagValue<'static>>);

    impl Rec {
        fn new(pairs: &[([u8; 2], TagValue<'static>)]) -> Self {
            Rec(pairs.iter().cloned().collect())
        }
    }

    impl TagSource for Rec {
        fn tag_value(&self, tag: [u8; 2]) -> anyhow::Result<TagValue<'_>> {
            Ok(self.0.get(&tag).cloned().unwrap_or(TagValue::Missing))
        }
    }

    fn s(v: &'static str) -> TagValue<'static> {
        TagValue::Str(Cow::Borrowed(v.as_bytes()))
    }

    fn keeps(expr: &str, rec: &Rec) -> bool {
        TagFilter::parse(expr).unwrap().keeps(rec).unwrap()
    }

    fn ont() -> Rec {
        Rec::new(&[
            (*b"er", s("signal_positive")),
            (*b"dx", TagValue::Num(1.0)),
            (*b"qs", TagValue::Num(12.5)),
            (*b"ch", TagValue::Num(1200.0)),
            (*b"st", s("2026-09-02T10:00:00Z")),
            (*b"mv", TagValue::Array),
        ])
    }

    #[test]
    fn every_documented_example_parses() {
        for e in [
            r#"[er]!="data_service_unblock_mux_change""#,
            r#"[er]=="signal_positive" || [er]=="signal_negative""#,
            "[dx]==1",
            "[dx]!=-1",
            "[qs]>=10 && [ch]<1500",
            "[rq]>=0.99 && [np]>=3",
            r#"[BC]=="barcode03" || ![BC]"#,
            "exists([fi])",
            "![MM]",
            "[dx]==1 || ([dx]==0 && [qs]>=15)",
            r#"[st]>="2026-09-01T00:00:00""#,
            "[qs]>=1e1",
        ] {
            TagFilter::parse(e).unwrap_or_else(|err| panic!("{e}: {err}"));
        }
    }

    #[test]
    fn comparisons_and_boolean_operators_follow_precedence() {
        let r = ont();
        assert!(keeps("[dx]==1", &r));
        assert!(!keeps("[dx]!=1", &r));
        assert!(keeps("[qs]>=10 && [ch]<1500", &r));
        assert!(!keeps("[qs]>=10 && [ch]<1000", &r));
        assert!(keeps("[qs]<10 || [ch]<1500", &r));
        assert!(
            keeps("[dx]==0 || [dx]==1 && [qs]>12", &r),
            "&& binds tighter than ||"
        );
        assert!(!keeps("([dx]==0 || [dx]==1) && [qs]>13", &r));
        assert!(keeps("!([dx]==0)", &r));
        assert!(keeps("!!([dx]==1)", &r));
        assert!(keeps(r#"[er]=="signal_positive""#, &r));
        assert!(keeps(r#"[st]>="2026-09-01T00:00:00""#, &r));
        assert!(!keeps(r#"[st]<"2026-09-01""#, &r));
    }

    #[test]
    fn presence_tests_cover_scalars_and_arrays() {
        let r = ont();
        assert!(keeps("[mv]", &r));
        assert!(keeps("exists([mv])", &r));
        assert!(!keeps("[BC]", &r));
        assert!(keeps("![BC]", &r));
        assert!(keeps(r#"[BC]=="barcode03" || ![BC]"#, &r));
    }

    /// A missing tag makes every comparison false, `!=` included, and only a
    /// presence test or a negated comparison accepts absence.
    #[test]
    fn missing_tags_fail_every_comparison() {
        let r = ont();
        assert!(!keeps("[BC]==\"x\"", &r));
        assert!(!keeps("[BC]!=\"x\"", &r));
        assert!(!keeps("[np]>=3", &r));
        assert!(!keeps("[np]<3", &r));
        assert!(keeps("!([np]>=3)", &r));
    }

    #[test]
    fn type_mismatches_are_errors_that_name_the_tag() {
        let r = ont();
        let f = TagFilter::parse("[qs]==\"10\"").unwrap();
        let err = f.keeps(&r).unwrap_err().to_string();
        assert!(err.contains("[qs]") && err.contains("a number"), "{err}");
        let f = TagFilter::parse("[mv]==1").unwrap();
        let err = f.keeps(&r).unwrap_err().to_string();
        assert!(err.contains("[mv]") && err.contains("an array"), "{err}");
        let f = TagFilter::parse("[er]<1").unwrap();
        assert!(f.keeps(&r).is_err());
    }

    #[test]
    fn integer_widths_and_floats_compare_alike() {
        let a = Rec::new(&[(*b"qs", TagValue::Num(10.0))]);
        let b = Rec::new(&[(*b"qs", TagValue::Num(10.0))]);
        assert!(keeps("[qs]>=10", &a) && keeps("[qs]>=10", &b));
        assert!(keeps("[qs]==10.0", &a));
    }

    #[test]
    fn parse_errors_name_the_problem() {
        for (expr, needle) in [
            ("[qs]", ""),
            ("10", "a literal is not a condition"),
            ("\"x\"", "a literal is not a condition"),
            ("[qs] >= ", "ends early"),
            ("[qs] >= 10 )", "unexpected )"),
            ("[q]==1", "invalid tag"),
            ("qs>=10", "unsupported word"),
            ("[qs]=~\"x\"", "unexpected character"),
            ("[er]==\"abc", "unterminated string"),
            ("[er]==\"a\\nb\"", "invalid escape"),
            ("length([qs])", "unsupported word"),
            ("10 == 10", "needs a tag"),
            ("![qs] == 1", ""),
        ] {
            let res = TagFilter::parse(expr);
            if needle.is_empty() {
                assert!(res.is_ok(), "{expr}: {res:?}");
            } else {
                let err = res.unwrap_err().to_string();
                assert!(err.contains(needle), "{expr}: {err}");
                assert!(err.contains("--tag-filter"), "{err}");
            }
        }
    }

    #[test]
    fn several_filters_all_have_to_hold() {
        let r = ont();
        let f = TagFilters::parse(&["[qs]>=10".to_string(), "[dx]==1".to_string()]).unwrap();
        assert!(f.keeps(&r).unwrap());
        let f = TagFilters::parse(&["[qs]>=10".to_string(), "[dx]==0".to_string()]).unwrap();
        assert!(!f.keeps(&r).unwrap());
        assert!(TagFilters::parse(&[]).unwrap().keeps(&r).unwrap());
        assert_eq!(
            TagFilters::parse(&["[qs]>=10".to_string()])
                .unwrap()
                .texts()
                .collect::<Vec<_>>(),
            ["[qs]>=10"]
        );
    }

    #[test]
    fn record_buf_data_is_a_tag_source() {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::Data;
        use noodles_sam::alignment::record_buf::data::field::Value;
        let mut data = Data::default();
        data.insert(Tag::new(b'd', b'x'), Value::Int8(1));
        data.insert(Tag::new(b'q', b's'), Value::Float(9.5));
        data.insert(
            Tag::new(b'e', b'r'),
            Value::String(b"signal_positive".as_slice().into()),
        );
        let f = TagFilters::parse(&["[dx]==1 && [qs]<10 && [er]!=\"x\"".to_string()]).unwrap();
        assert!(f.keeps(&data).unwrap());
    }
}
