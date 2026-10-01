//! Keyed output routing: the per-key parts of a rendered batch and the output
//! sinks they are written to, one per key.

/// Identifies one output of a run. A single-output run writes key 0 only.
pub(crate) type KeyId = u32;

/// One batch's rendered output, split by key. Each part keeps its items in
/// the order they were pushed. Key 0 is held inline, so a batch that renders
/// to key 0 alone allocates exactly one buffer.
#[derive(Clone)]
pub(crate) struct Parts<T> {
    /// The part for key 0.
    first: Vec<T>,
    /// The parts for keys 1 and up, at index `key - 1`.
    rest: Vec<Vec<T>>,
}

impl<T> Default for Parts<T> {
    fn default() -> Self {
        Self {
            first: Vec::new(),
            rest: Vec::new(),
        }
    }
}

impl<T> Parts<T> {
    /// Returns the part for `key`, creating it empty on first use.
    pub(crate) fn part(&mut self, key: KeyId) -> &mut Vec<T> {
        let Some(idx) = (key as usize).checked_sub(1) else {
            return &mut self.first;
        };
        if idx >= self.rest.len() {
            self.rest.resize_with(idx + 1, Vec::new);
        }
        &mut self.rest[idx]
    }

    /// Whether every part other than key 0's is empty.
    pub(crate) fn only_first(&self) -> bool {
        self.rest.iter().all(Vec::is_empty)
    }

    /// Yields every non-empty part with its key, in ascending key order.
    pub(crate) fn into_nonempty(self) -> impl Iterator<Item = (KeyId, Vec<T>)> {
        std::iter::once(self.first)
            .chain(self.rest)
            .enumerate()
            .filter(|(_, part)| !part.is_empty())
            .map(|(key, part)| (key as KeyId, part))
    }

    /// Hands every non-empty part to `f` with its key, in ascending key
    /// order, and empties it, keeping its allocation for reuse. Stops at the
    /// first error, leaving the remaining parts as they are.
    pub(crate) fn drain_each(
        &mut self,
        mut f: impl FnMut(KeyId, &[T]) -> anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        let parts = std::iter::once(&mut self.first).chain(self.rest.iter_mut());
        for (key, part) in parts.enumerate() {
            if !part.is_empty() {
                f(key as KeyId, part)?;
                part.clear();
            }
        }
        Ok(())
    }
}

/// Opens the sink for one key.
pub(crate) type Opener<S> = Box<dyn FnMut(KeyId) -> anyhow::Result<S> + Send>;

/// The output sinks of a run, indexed by key. A sink is opened on the first
/// `get` for its key and kept until `into_sinks` hands it back for
/// finalizing.
pub(crate) struct KeyedSinks<S> {
    /// The sink for each key, at index `key`; `None` until opened.
    sinks: Vec<Option<S>>,
    /// Opens a sink for a key not yet open; `None` for a fixed set.
    open: Option<Opener<S>>,
}

impl<S> KeyedSinks<S> {
    /// A single-output set: `sink` is key 0, and no other key can be opened.
    pub(crate) fn single(sink: S) -> Self {
        Self {
            sinks: vec![Some(sink)],
            open: None,
        }
    }

    /// An empty set that opens each key's sink with `open` on first use.
    pub(crate) fn with_opener(open: Opener<S>) -> Self {
        Self {
            sinks: Vec::new(),
            open: Some(open),
        }
    }

    /// Returns the sink for `key`, opening it on first use. An opener error
    /// leaves the key unopened; a key outside a single-output set is an error.
    pub(crate) fn get(&mut self, key: KeyId) -> anyhow::Result<&mut S> {
        let idx = key as usize;
        if idx >= self.sinks.len() {
            self.sinks.resize_with(idx + 1, || None);
        }
        let slot = &mut self.sinks[idx];
        if slot.is_none() {
            let Some(open) = self.open.as_mut() else {
                anyhow::bail!("no output is configured for key {key}");
            };
            *slot = Some(open(key)?);
        }
        Ok(slot.as_mut().expect("The slot was filled above"))
    }

    /// Returns every opened sink with its key, in ascending key order.
    pub(crate) fn into_sinks(self) -> Vec<(KeyId, S)> {
        self.sinks
            .into_iter()
            .enumerate()
            .filter_map(|(key, sink)| sink.map(|sink| (key as KeyId, sink)))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    use super::*;

    #[test]
    fn parts_keep_per_key_order() {
        let mut parts: Parts<char> = Parts::default();
        parts.part(2).push('a');
        parts.part(0).push('b');
        parts.part(2).push('c');
        let out: Vec<(KeyId, Vec<char>)> = parts.into_nonempty().collect();
        assert_eq!(out, [(0, vec!['b']), (2, vec!['a', 'c'])]);
    }

    #[test]
    fn drain_each_visits_in_key_order_and_empties() {
        let mut parts: Parts<u8> = Parts::default();
        parts.part(2).push(20);
        parts.part(0).extend([1, 2]);
        let mut seen = Vec::new();
        parts
            .drain_each(|key, part| {
                seen.push((key, part.to_vec()));
                Ok(())
            })
            .unwrap();
        assert_eq!(seen, [(0, vec![1, 2]), (2, vec![20])]);
        assert_eq!(parts.into_nonempty().count(), 0);
    }

    #[test]
    fn parts_skip_empty_keys() {
        let mut parts: Parts<u8> = Parts::default();
        parts.part(3);
        assert_eq!(parts.into_nonempty().count(), 0);
    }

    #[test]
    fn keyed_sinks_open_lazily_once() {
        let opened = Arc::new(AtomicUsize::new(0));
        let counter = Arc::clone(&opened);
        let mut sinks: KeyedSinks<Vec<KeyId>> = KeyedSinks::with_opener(Box::new(move |key| {
            counter.fetch_add(1, Ordering::Relaxed);
            Ok(vec![key])
        }));
        sinks.get(3).unwrap().push(30);
        sinks.get(3).unwrap().push(31);
        sinks.get(1).unwrap().push(10);
        assert_eq!(opened.load(Ordering::Relaxed), 2);
        assert_eq!(sinks.into_sinks(), [(1, vec![1, 10]), (3, vec![3, 30, 31])]);
    }

    #[test]
    fn single_sink_holds_key_zero_only() {
        let mut sinks = KeyedSinks::single(vec![7u8]);
        sinks.get(0).unwrap().push(8);
        assert!(sinks.get(1).is_err());
        assert_eq!(sinks.into_sinks(), [(0, vec![7, 8])]);
    }

    #[test]
    fn opener_error_leaves_the_key_unopened() {
        let mut sinks: KeyedSinks<Vec<u8>> =
            KeyedSinks::with_opener(Box::new(|key| anyhow::bail!("cannot open key {key}")));
        assert_eq!(sinks.get(2).unwrap_err().to_string(), "cannot open key 2");
        assert!(sinks.into_sinks().is_empty());
    }
}
