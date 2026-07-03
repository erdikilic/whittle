/// Format-neutral read carrier. `qual` holds raw Phred scores (0-based), not ASCII.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadRecord {
    pub name: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}
