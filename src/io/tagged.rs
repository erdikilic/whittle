//! FASTQ headers that carry SAM auxiliary fields.
//!
//! `samtools fastq -T` and whittle's own BAM-to-FASTQ output append aux tags to
//! the header line, tab-delimited and spelled as in SAM text
//! (`@name\tMM:Z:C+m,3;\tML:B:C,200`). A read from such a file is decoded into a
//! `RecordBuf` and takes the BAM-to-FASTQ path, so its tags are rewritten per
//! output segment exactly as a uBAM record's are.

use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::Data;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;

use crate::record::ReadRecord;

/// Returns true when `head` (a FASTQ header line without the leading `@`)
/// carries SAM aux fields: a tab followed by a two-character tag, a colon, a
/// type character and a colon.
pub fn has_aux_tags(head: &[u8]) -> bool {
    let Some(tab) = head.iter().position(|&b| b == b'\t') else {
        return false;
    };
    let f = &head[tab + 1..];
    f.len() >= 5
        && f[0].is_ascii_alphabetic()
        && f[1].is_ascii_alphanumeric()
        && f[2] == b':'
        && matches!(f[3], b'A' | b'i' | b'f' | b'Z' | b'H' | b'B')
        && f[4] == b':'
}

/// Decodes a tagged FASTQ read into an unmapped `RecordBuf`: the name is the
/// first whitespace-delimited token, and every tab-delimited field is an
/// aux tag. A field that does not parse names the read and the field.
pub fn record_from_tagged(rec: ReadRecord) -> anyhow::Result<RecordBuf> {
    let (name, fields) = match rec.name.iter().position(|&b| b == b'\t') {
        Some(tab) => (&rec.name[..tab], &rec.name[tab + 1..]),
        None => (rec.name.as_slice(), &[][..]),
    };
    let data = parse_aux_fields(fields)
        .map_err(|e| anyhow::anyhow!("read {}: header {e}", String::from_utf8_lossy(name)))?;
    let mut out = RecordBuf::default();
    *out.flags_mut() = Flags::UNMAPPED;
    let id_end = name
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(name.len());
    *out.name_mut() = Some(name[..id_end].to_vec().into());
    *out.sequence_mut() = rec.seq.into();
    *out.quality_scores_mut() = rec.qual.into();
    *out.data_mut() = data;
    Ok(out)
}

/// Parses tab-delimited SAM aux fields (`TAG:TYPE:VALUE`) into `Data`. Integer
/// values take the narrowest BAM subtype that holds them, as a BAM reader
/// would present them. Empty input yields empty data.
pub fn parse_aux_fields(fields: &[u8]) -> anyhow::Result<Data> {
    let mut data = Data::default();
    if fields.is_empty() {
        return Ok(data);
    }
    for field in fields.split(|&b| b == b'\t') {
        let (tag, value) = parse_aux_field(field)
            .map_err(|e| anyhow::anyhow!("field {:?}: {e}", String::from_utf8_lossy(field)))?;
        if data.insert(tag, value).is_some() {
            anyhow::bail!("field {:?}: duplicate tag", String::from_utf8_lossy(field));
        }
    }
    Ok(data)
}

/// Parses one `TAG:TYPE:VALUE` field.
fn parse_aux_field(field: &[u8]) -> anyhow::Result<(Tag, Value)> {
    if field.len() < 5 || field[2] != b':' || field[4] != b':' {
        anyhow::bail!("expected TAG:TYPE:VALUE");
    }
    if !(field[0].is_ascii_alphabetic() && field[1].is_ascii_alphanumeric()) {
        anyhow::bail!("tag must be a letter followed by a letter or digit");
    }
    let tag = Tag::from([field[0], field[1]]);
    let raw = &field[5..];
    let value = match field[3] {
        b'A' => match raw {
            [c] if (b'!'..=b'~').contains(c) => Value::Character(*c),
            _ => anyhow::bail!("type A takes one printable character"),
        },
        b'i' => {
            let n: i64 =
                parse_int(raw, true).ok_or_else(|| anyhow::anyhow!("type i takes an integer"))?;
            Value::try_from(n).map_err(|_| anyhow::anyhow!("integer {n} exceeds 32 bits"))?
        },
        b'f' => {
            Value::Float(parse_float(raw).ok_or_else(|| anyhow::anyhow!("type f takes a number"))?)
        },
        b'Z' => {
            if !raw.iter().all(|b| (b' '..=b'~').contains(b)) {
                anyhow::bail!("type Z takes printable text");
            }
            Value::String(raw.to_vec().into())
        },
        b'H' => {
            if !raw.len().is_multiple_of(2) || !raw.iter().all(u8::is_ascii_hexdigit) {
                anyhow::bail!("type H takes an even number of hex digits");
            }
            Value::Hex(raw.to_vec().into())
        },
        b'B' => Value::Array(parse_array(raw)?),
        other => anyhow::bail!("unknown type {:?}", char::from(other)),
    };
    Ok((tag, value))
}

/// Parses the value of a `B` field: a subtype character, then zero or more
/// comma-separated numbers.
fn parse_array(raw: &[u8]) -> anyhow::Result<Array> {
    let Some((&subtype, rest)) = raw.split_first() else {
        anyhow::bail!("type B takes a subtype");
    };
    let items = match rest {
        [] => &[][..],
        [b',', items @ ..] => items,
        _ => anyhow::bail!("type B takes a subtype followed by comma-separated values"),
    };
    /// Parses the comma-separated elements of an integer array; `signed`
    /// admits a leading `-`.
    fn ints<T: TryFrom<i64>>(items: &[u8], signed: bool) -> anyhow::Result<Vec<T>> {
        if items.is_empty() {
            return Ok(Vec::new());
        }
        let mut out = Vec::with_capacity(items.iter().filter(|&&b| b == b',').count() + 1);
        let mut rest = items;
        loop {
            // A one-byte element, the form of every move-table entry, is split
            // off without a search for its comma.
            let (item, tail) = match rest {
                [b, b',', tail @ ..] if *b != b',' => (&rest[..1], Some(tail)),
                _ => match rest.iter().position(|&b| b == b',') {
                    Some(i) => (&rest[..i], Some(&rest[i + 1..])),
                    None => (rest, None),
                },
            };
            out.push(parse_int(item, signed).ok_or_else(|| {
                anyhow::anyhow!(
                    "array value {:?} does not fit the subtype",
                    String::from_utf8_lossy(item)
                )
            })?);
            match tail {
                Some(tail) => rest = tail,
                None => return Ok(out),
            }
        }
    }
    Ok(match subtype {
        b'c' => Array::Int8(ints(items, true)?),
        b'C' => Array::UInt8(ints(items, false)?),
        b's' => Array::Int16(ints(items, true)?),
        b'S' => Array::UInt16(ints(items, false)?),
        b'i' => Array::Int32(ints(items, true)?),
        b'I' => Array::UInt32(ints(items, false)?),
        b'f' => {
            if items.is_empty() {
                Array::Float(Vec::new())
            } else {
                Array::Float(
                    items
                        .split(|&b| b == b',')
                        .map(|s| {
                            parse_float(s).ok_or_else(|| {
                                anyhow::anyhow!(
                                    "array value {:?} is not a number",
                                    String::from_utf8_lossy(s)
                                )
                            })
                        })
                        .collect::<anyhow::Result<_>>()?,
                )
            }
        },
        other => anyhow::bail!("unknown array subtype {:?}", char::from(other)),
    })
}

/// Parses a decimal integer that fits `T`, accepting what `str::parse` accepts:
/// an optional `+`, or `-` when `signed`, then one or more ASCII digits.
fn parse_int<T: TryFrom<i64>>(raw: &[u8], signed: bool) -> Option<T> {
    let (negative, digits) = match raw {
        [b'+', rest @ ..] => (false, rest),
        [b'-', rest @ ..] if signed => (true, rest),
        _ => (false, raw),
    };
    if digits.is_empty() {
        return None;
    }
    let mut magnitude = 0u64;
    for &b in digits {
        let digit = b.wrapping_sub(b'0');
        if digit > 9 {
            return None;
        }
        magnitude = magnitude.checked_mul(10)?.checked_add(u64::from(digit))?;
    }
    let value = if negative {
        0i64.checked_sub_unsigned(magnitude)?
    } else {
        i64::try_from(magnitude).ok()?
    };
    T::try_from(value).ok()
}

/// Parses a SAM float: a decimal or exponent form, `nan`, `inf`, `-inf`.
fn parse_float(raw: &[u8]) -> Option<f32> {
    std::str::from_utf8(raw).ok()?.parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tagged(name: &str, fields: &str) -> ReadRecord {
        ReadRecord {
            name: format!("{name}\t{fields}").into_bytes(),
            seq: b"ACGT".to_vec(),
            qual: vec![30; 4],
        }
    }

    #[test]
    fn detects_a_sam_field_after_the_first_tab() {
        assert!(has_aux_tags(b"r1\tMM:Z:C+m,1;"));
        assert!(has_aux_tags(b"r1\tqs:f:12.5\tML:B:C,1"));
        assert!(has_aux_tags(b"r1 desc\tRG:Z:x"));
        assert!(!has_aux_tags(b"r1"));
        assert!(!has_aux_tags(b"r1 runid=abc sampleid=x"));
        assert!(!has_aux_tags(b"r1\tcomment"));
        assert!(!has_aux_tags(b"r1\tMM:X:1"));
        assert!(!has_aux_tags(b"r1\t"));
    }

    #[test]
    fn decodes_every_field_type() {
        let rec = tagged(
            "r1",
            "MM:Z:C+m,1;\tML:B:C,200,7\tMN:i:4\tqs:f:12.5\tXA:A:x\tXH:H:1A2B\tXS:B:s,-3,4\tXF:B:f,0.5,1e3\tXE:B:I",
        );
        let out = record_from_tagged(rec).unwrap();
        assert_eq!(out.name().map(|n| n.as_ref()), Some(&b"r1"[..]));
        assert_eq!(out.sequence().as_ref(), b"ACGT");
        assert_eq!(out.quality_scores().as_ref(), &[30; 4]);
        assert!(out.flags().is_unmapped());
        let d = out.data();
        assert_eq!(
            d.get(&Tag::BASE_MODIFICATIONS),
            Some(&Value::String(b"C+m,1;".to_vec().into()))
        );
        assert_eq!(
            d.get(&Tag::BASE_MODIFICATION_PROBABILITIES),
            Some(&Value::Array(Array::UInt8(vec![200, 7])))
        );
        assert_eq!(d.get(&Tag::from(*b"MN")), Some(&Value::UInt8(4)));
        assert_eq!(d.get(&Tag::from(*b"qs")), Some(&Value::Float(12.5)));
        assert_eq!(d.get(&Tag::from(*b"XA")), Some(&Value::Character(b'x')));
        assert_eq!(
            d.get(&Tag::from(*b"XH")),
            Some(&Value::Hex(b"1A2B".to_vec().into()))
        );
        assert_eq!(
            d.get(&Tag::from(*b"XS")),
            Some(&Value::Array(Array::Int16(vec![-3, 4])))
        );
        assert_eq!(
            d.get(&Tag::from(*b"XF")),
            Some(&Value::Array(Array::Float(vec![0.5, 1000.0])))
        );
        assert_eq!(
            d.get(&Tag::from(*b"XE")),
            Some(&Value::Array(Array::UInt32(Vec::new())))
        );
    }

    #[test]
    fn integers_take_the_narrowest_subtype() {
        let d = parse_aux_fields(b"a1:i:-5\ta2:i:300\ta3:i:-40000\ta4:i:70000").unwrap();
        assert_eq!(d.get(&Tag::from(*b"a1")), Some(&Value::Int8(-5)));
        assert_eq!(d.get(&Tag::from(*b"a2")), Some(&Value::UInt16(300)));
        assert_eq!(d.get(&Tag::from(*b"a3")), Some(&Value::Int32(-40000)));
        assert_eq!(d.get(&Tag::from(*b"a4")), Some(&Value::UInt32(70000)));
    }

    #[test]
    fn a_header_without_tags_decodes_with_empty_data() {
        let out = record_from_tagged(ReadRecord {
            name: b"r1".to_vec(),
            seq: b"AC".to_vec(),
            qual: vec![1, 2],
        })
        .unwrap();
        assert!(out.data().is_empty());
        assert_eq!(out.name().map(|n| n.as_ref()), Some(&b"r1"[..]));
    }

    /// The byte-level integer parser accepts and rejects exactly what
    /// `str::parse` does for every integer width the arrays and `i` fields use.
    #[test]
    fn integer_parser_matches_str_parse() {
        fn check<T: TryFrom<i64> + std::str::FromStr + PartialEq + std::fmt::Debug>(signed: bool) {
            for text in [
                "",
                "+",
                "-",
                "0",
                "-0",
                "+0",
                "7",
                "+7",
                "-7",
                "007",
                "-007",
                "127",
                "128",
                "-128",
                "-129",
                "255",
                "256",
                "32767",
                "32768",
                "-32768",
                "-32769",
                "65535",
                "65536",
                "2147483647",
                "2147483648",
                "-2147483648",
                "-2147483649",
                "4294967295",
                "4294967296",
                "9223372036854775807",
                "9223372036854775808",
                "-9223372036854775808",
                "-9223372036854775809",
                "18446744073709551615",
                "18446744073709551616",
                "000000000000000000000000000001",
                "1 ",
                " 1",
                "1a",
                "0x10",
                "--1",
                "+-1",
                "-+1",
                "1.0",
                "1e3",
                "\u{661}",
            ] {
                assert_eq!(
                    parse_int::<T>(text.as_bytes(), signed),
                    text.parse::<T>().ok(),
                    "{text:?} as {}",
                    std::any::type_name::<T>()
                );
            }
        }
        check::<i8>(true);
        check::<u8>(false);
        check::<i16>(true);
        check::<u16>(false);
        check::<i32>(true);
        check::<u32>(false);
        check::<i64>(true);
    }

    #[test]
    fn malformed_fields_name_the_read_and_the_field() {
        for (fields, msg) in [
            ("MM:Z", "expected TAG:TYPE:VALUE"),
            ("1M:i:3", "tag must be a letter"),
            ("XX:i:abc", "type i takes an integer"),
            ("XX:i:99999999999", "exceeds 32 bits"),
            ("XX:Q:1", "unknown type"),
            ("XX:B:C,300", "does not fit the subtype"),
            ("XX:B:C,1,-0", "does not fit the subtype"),
            ("XX:B:c,1,", "does not fit the subtype"),
            ("XX:B:c,1,,2", "does not fit the subtype"),
            ("XX:B:c,,", "does not fit the subtype"),
            ("XX:i:+", "type i takes an integer"),
            ("XX:B:x,1", "unknown array subtype"),
            ("XX:B:C;1", "comma-separated"),
            ("XX:H:ABC", "even number of hex digits"),
            ("XX:A:ab", "one printable character"),
            ("XX:i:1\tXX:i:2", "duplicate tag"),
        ] {
            let err = record_from_tagged(tagged("r7", fields))
                .unwrap_err()
                .to_string();
            assert!(err.starts_with("read r7: header field"), "{err}");
            assert!(err.contains(msg), "{fields}: {err}");
        }
    }
}
