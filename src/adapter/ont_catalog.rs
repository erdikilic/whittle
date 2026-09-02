//! The built-in ONT adapter, primer, and barcode catalog.
//!
//! Assembled for whittle on 2026-07-06 from primary sources, cross-verified:
//! dorado `adapter_primer_kits.cpp` and `utils/barcode_kits.cpp` (kit-14
//! authoritative), Porechop `porechop/adapters.py` (legacy plus the 96
//! barcodes), and the qcat kit YAMLs. The sequences are ONT-published facts;
//! the catalog itself is whittle's own compilation.
//!
//! `End::Five` means the sequence is expected at the read start, `End::Three` at
//! the read end, `End::Both` at either. Every entry is also searched
//! reverse-complemented, so opposite-strand reads are covered regardless.
//!
//! Barcode dedup rule: one barcode number is a single shared 24 bp oligo across
//! the native, PCR, and rapid kits, which differ only in flank and orientation
//! (native uses the reverse complement). Only the canonical forward oligo is
//! stored, and reverse-complement search covers the rest. Trimming through the
//! flanks removes the barcode whatever its number.
//!
//! Each entry's comment names the kits it belongs to and the upstream source it
//! was taken from.

use super::End;

/// One catalog entry: display name, the end it is expected at, and the sequence.
/// Sequences are uppercase nucleotide codes of at least `MIN_PATTERN_LEN`
/// bases, which `preset::tests::entries_are_valid_nucleotide_sequences` and
/// `preset::tests::entries_meet_the_minimum_pattern_length` enforce.
pub(super) type Entry = (&'static str, End, &'static [u8]);

#[rustfmt::skip]
pub(super) const CATALOG: &[Entry] = &[

    // Ligation Y-adapter, two chemistry generations, both kept.
    // SQK-LSK114 and all *114 (PCS/RAD/ULK/NBD/PCB/RBK/16S)
    // [dorado:adapter_primer_kits.cpp(LSK110)]
    ("LSK114_front", End::Five, b"CCTGTACTTCGTTCAGTTACGTATTGC"),
    // SQK-LSK114 and all *114 [dorado:adapter_primer_kits.cpp(LSK110)]
    ("LSK114_rear", End::Three, b"AGCAATACGTAACTGAAC"),
    // SQK-LSK108/109, NSK007 (kit 9/10) [porechop:adapters.py(SQK-NSK007_Y_Top)]
    ("LSK109_front", End::Five, b"AATGTACTTCGTTCAGTTACGTATTGCT"),
    // SQK-LSK108/109, NSK007 (kit 9/10) [porechop:adapters.py(SQK-NSK007_Y_Bottom)]
    ("LSK109_rear", End::Three, b"GCAATACGTAACTGAACGAAGT"),

    // Rapid adapter, which doubles as the 3' flank of the rapid-barcoding kits.
    // SQK-RAD004/RAD114, ULK114, RBK* [dorado:adapter_primer_kits.cpp(RAD) /
    // porechop(Rapid_adapter)]
    ("RAD", End::Both, b"GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"),

    // Direct RNA.
    // SQK-RNA004, RNA004-XL, DRB004 [dorado:adapter_primer_kits.cpp(RNA004)]
    ("RNA004_rear", End::Three, b"GGTTGTTTCTGTTGGTGCTG"),

    // Primers: PCR, cDNA, and 10X.
    // SQK-LSK114 cDNA [dorado:adapter_primer_kits.cpp(cDNA)]
    ("cDNA_front", End::Five, b"TTTCTGTTGGTGCTGATATTGCTGGG"),
    // SQK-LSK114 cDNA [dorado:adapter_primer_kits.cpp(cDNA)]
    ("cDNA_rear", End::Three, b"ACTTGCCTGTCGCTCTATCTTCTTT"),
    // SQK-PCS114, PCB114 [dorado:adapter_primer_kits.cpp(PCS110)]
    ("PCS110_front", End::Five, b"TTTCTGTTGGTGCTGATATTGCTTT"),
    // SQK-PCS114, PCB114 [dorado:adapter_primer_kits.cpp(PCS110)]
    ("PCS110_rear", End::Three, b"ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAGTTTT"),
    // SQK-LSK114 (10X) [dorado:adapter_primer_kits.cpp(GEN10X)]
    ("GEN10X_front", End::Five, b"CTACACGACGCTCTTCCGATCT"),
    // SQK-LSK114 (10X) [dorado:adapter_primer_kits.cpp(GEN10X)]
    ("GEN10X_rear", End::Three, b"GTACTCTGCGTTGATACCACTGCTT"),
    // PCR/cDNA (legacy) [porechop:adapters.py(PCR_1)]
    ("PCR1_front", End::Five, b"ACTTGCCTGTCGCTCTATCTTC"),
    // PCR/cDNA (legacy) [porechop:adapters.py(PCR_1)]
    ("PCR1_rear", End::Three, b"GAAGATAGAGCGACAGGCAAGT"),
    // PCR/cDNA (legacy) [porechop:adapters.py(PCR_2)]
    ("PCR2_front", End::Five, b"TTTCTGTTGGTGCTGATATTGC"),
    // PCR/cDNA (legacy) [porechop:adapters.py(PCR_2)]
    ("PCR2_rear", End::Three, b"GCAATATCAGCACCAACAGAAA"),
    // cDNA SSP (legacy) [porechop:adapters.py(cDNA_SSP)]
    ("cDNA_SSP", End::Five, b"TTTCTGTTGGTGCTGATATTGCTGCCATTACGGCCGGG"),

    // Barcode flanks (dorado kit-14 constants). Trimming through a flank removes the
    // barcode regardless of its number. Flanks shorter than `MIN_PATTERN_LEN`
    // (dorado's NB_1st_REAR, BC_1st_FRONT, RBK_FRONT and RLB_FRONT, 7 to 8 bp)
    // are omitted: a pattern that short is never searched standalone.
    // native barcoding (NBD*) [dorado:barcode_kits.cpp(NB_1st_FRONT)]
    ("NB_front", End::Five, b"ATTGCTAAGGTTAA"),
    // PCR barcoding (PBC/BC*) [dorado:barcode_kits.cpp(BC_1st_REAR)]
    ("PBC_rear", End::Three, b"TTAACCTTTCTGTTGGTGCTGATATTGC"),
    // rapid barcoding v4 / kit14 (RBK*) [dorado:barcode_kits.cpp(RBK4_FRONT)]
    ("RBK4_front", End::Five, b"GCTTGGGTGTTTAACC"),
    // rapid lig barcoding (RLB) [dorado:barcode_kits.cpp(RLB_REAR)]
    ("RLB_rear", End::Three, b"CGTTTTTCGTGCGCCGCTTC"),
    // 16S barcoding (RAB/16S) [dorado:barcode_kits.cpp(RAB_1st_REAR)]
    ("RAB_16S_rear1", End::Three, b"AGAGTTTGATCATGGCTCAG"),
    // 16S barcoding (RAB/16S) [dorado:barcode_kits.cpp(RAB_2nd_REAR)]
    ("RAB_16S_rear2", End::Three, b"CGGTTACCTTGTTACGACTT"),
    // ligation/PCR barcoding (LWB) [dorado:barcode_kits.cpp(LWB_1st_REAR)]
    ("LWB_rear1", End::Three, b"ACTTGCCTGTCGCTCTATCTTC"),
    // MAB [dorado:barcode_kits.cpp(MAB_REAR)]
    ("MAB_rear", End::Three, b"CCATATCCGTGTCGCCCTT"),

    // The 96 canonical forward barcode oligos, from Porechop `adapters.py` and confirmed
    // present in dorado. Searched at both ends, in both orientations.
    ("BC01", End::Both, b"AAGAAAGTTGTCGGTGTCTTTGTG"),
    ("BC02", End::Both, b"TCGATTCCGTTTGTAGTCGTCTGT"),
    ("BC03", End::Both, b"GAGTCTTGTGTCCCAGTTACCAGG"),
    ("BC04", End::Both, b"TTCGGATTCTATCGTGTTTCCCTA"),
    ("BC05", End::Both, b"CTTGTCCAGGGTTTGTGTAACCTT"),
    ("BC06", End::Both, b"TTCTCGCAAAGGCAGAAAGTAGTC"),
    ("BC07", End::Both, b"GTGTTACCGTGGGAATGAATCCTT"),
    ("BC08", End::Both, b"TTCAGGGAACAAACCAAGTTACGT"),
    ("BC09", End::Both, b"AACTAGGCACAGCGAGTCTTGGTT"),
    ("BC10", End::Both, b"AAGCGTTGAAACCTTTGTCCTCTC"),
    ("BC11", End::Both, b"GTTTCATCTATCGGAGGGAATGGA"),
    ("BC12", End::Both, b"CAGGTAGAAAGAAGCAGAATCGGA"),
    ("BC13", End::Both, b"AGAACGACTTCCATACTCGTGTGA"),
    ("BC14", End::Both, b"AACGAGTCTCTTGGGACCCATAGA"),
    ("BC15", End::Both, b"AGGTCTACCTCGCTAACACCACTG"),
    ("BC16", End::Both, b"CGTCAACTGACAGTGGTTCGTACT"),
    ("BC17", End::Both, b"ACCCTCCAGGAAAGTACCTCTGAT"),
    ("BC18", End::Both, b"CCAAACCCAACAACCTAGATAGGC"),
    ("BC19", End::Both, b"GTTCCTCGTGCAGTGTCAAGAGAT"),
    ("BC20", End::Both, b"TTGCGTCCTGTTACGAGAACTCAT"),
    ("BC21", End::Both, b"GAGCCTCTCATTGTCCGTTCTCTA"),
    ("BC22", End::Both, b"ACCACTGCCATGTATCAAAGTACG"),
    ("BC23", End::Both, b"CTTACTACCCAGTGAACCTCCTCG"),
    ("BC24", End::Both, b"GCATAGTTCTGCATGATGGGTTAG"),
    ("BC25", End::Both, b"GTAAGTTGGGTATGCAACGCAATG"),
    ("BC26", End::Both, b"CATACAGCGACTACGCATTCTCAT"),
    ("BC27", End::Both, b"CGACGGTTAGATTCACCTCTTACA"),
    ("BC28", End::Both, b"TGAAACCTAAGAAGGCACCGTATC"),
    ("BC29", End::Both, b"CTAGACACCTTGGGTTGACAGACC"),
    ("BC30", End::Both, b"TCAGTGAGGATCTACTTCGACCCA"),
    ("BC31", End::Both, b"TGCGTACAGCAATCAGTTACATTG"),
    ("BC32", End::Both, b"CCAGTAGAAGTCCGACAACGTCAT"),
    ("BC33", End::Both, b"CAGACTTGGTACGGTTGGGTAACT"),
    ("BC34", End::Both, b"GGACGAAGAACTCAAGTCAAAGGC"),
    ("BC35", End::Both, b"CTACTTACGAAGCTGAGGGACTGC"),
    ("BC36", End::Both, b"ATGTCCCAGTTAGAGGAGGAAACA"),
    ("BC37", End::Both, b"GCTTGCGATTGATGCTTAGTATCA"),
    ("BC38", End::Both, b"ACCACAGGAGGACGATACAGAGAA"),
    ("BC39", End::Both, b"CCACAGTGTCAACTAGAGCCTCTC"),
    ("BC40", End::Both, b"TAGTTTGGATGACCAAGGATAGCC"),
    ("BC41", End::Both, b"GGAGTTCGTCCAGAGAAGTACACG"),
    ("BC42", End::Both, b"CTACGTGTAAGGCATACCTGCCAG"),
    ("BC43", End::Both, b"CTTTCGTTGTTGACTCGACGGTAG"),
    ("BC44", End::Both, b"AGTAGAAAGGGTTCCTTCCCACTC"),
    ("BC45", End::Both, b"GATCCAACAGAGATGCCTTCAGTG"),
    ("BC46", End::Both, b"GCTGTGTTCCACTTCATTCTCCTG"),
    ("BC47", End::Both, b"GTGCAACTTTCCCACAGGTAGTTC"),
    ("BC48", End::Both, b"CATCTGGAACGTGGTACACCTGTA"),
    ("BC49", End::Both, b"ACTGGTGCAGCTTTGAACATCTAG"),
    ("BC50", End::Both, b"ATGGACTTTGGTAACTTCCTGCGT"),
    ("BC51", End::Both, b"GTTGAATGAGCCTACTGGGTCCTC"),
    ("BC52", End::Both, b"TGAGAGACAAGATTGTTCGTGGAC"),
    ("BC53", End::Both, b"AGATTCAGACCGTCTCATGCAAAG"),
    ("BC54", End::Both, b"CAAGAGCTTTGACTAAGGAGCATG"),
    ("BC55", End::Both, b"TGGAAGATGAGACCCTGATCTACG"),
    ("BC56", End::Both, b"TCACTACTCAACAGGTGGCATGAA"),
    ("BC57", End::Both, b"GCTAGGTCAATCTCCTTCGGAAGT"),
    ("BC58", End::Both, b"CAGGTTACTCCTCCGTGAGTCTGA"),
    ("BC59", End::Both, b"TCAATCAAGAAGGGAAAGCAAGGT"),
    ("BC60", End::Both, b"CATGTTCAACCAAGGCTTCTATGG"),
    ("BC61", End::Both, b"AGAGGGTACTATGTGCCTCAGCAC"),
    ("BC62", End::Both, b"CACCCACACTTACTTCAGGACGTA"),
    ("BC63", End::Both, b"TTCTGAAGTTCCTGGGTCTTGAAC"),
    ("BC64", End::Both, b"GACAGACACCGTTCATCGACTTTC"),
    ("BC65", End::Both, b"TTCTCAGTCTTCCTCCAGACAAGG"),
    ("BC66", End::Both, b"CCGATCCTTGTGGCTTCTAACTTC"),
    ("BC67", End::Both, b"GTTTGTCATACTCGTGTGCTCACC"),
    ("BC68", End::Both, b"GAATCTAAGCAAACACGAAGGTGG"),
    ("BC69", End::Both, b"TACAGTCCGAGCCTCATGTGATCT"),
    ("BC70", End::Both, b"ACCGAGATCCTACGAATGGAGTGT"),
    ("BC71", End::Both, b"CCTGGGAGCATCAGGTAGTAACAG"),
    ("BC72", End::Both, b"TAGCTGACTGTCTTCCATACCGAC"),
    ("BC73", End::Both, b"AAGAAACAGGATGACAGAACCCTC"),
    ("BC74", End::Both, b"TACAAGCATCCCAACACTTCCACT"),
    ("BC75", End::Both, b"GACCATTGTGATGAACCCTGTTGT"),
    ("BC76", End::Both, b"ATGCTTGTTACATCAACCCTGGAC"),
    ("BC77", End::Both, b"CGACCTGTTTCTCAGGGATACAAC"),
    ("BC78", End::Both, b"AACAACCGAACCTTTGAATCAGAA"),
    ("BC79", End::Both, b"TCTCGGAGATAGTTCTCACTGCTG"),
    ("BC80", End::Both, b"CGGATGAACATAGGATAGCGATTC"),
    ("BC81", End::Both, b"CCTCATCTTGTGAAGTTGTTTCGG"),
    ("BC82", End::Both, b"ACGGTATGTCGAGTTCCAGGACTA"),
    ("BC83", End::Both, b"TGGCTTGATCTAGGTAAGGTCGAA"),
    ("BC84", End::Both, b"GTAGTGGACCTAGAACCTGTGCCA"),
    ("BC85", End::Both, b"AACGGAGGAGTTAGTTGGATGATC"),
    ("BC86", End::Both, b"AGGTGATCCCAACAAGCGTAAGTA"),
    ("BC87", End::Both, b"TACATGCTCCTGTTGTTAGGGAGG"),
    ("BC88", End::Both, b"TCTTCTACTACCGATCCGAAGCAG"),
    ("BC89", End::Both, b"ACAGCATCAATGTTTGGCTAGTTG"),
    ("BC90", End::Both, b"GATGTAGAGGGTACGGTTTGAGGC"),
    ("BC91", End::Both, b"GGCTCCATAGGAACTCACGCTACT"),
    ("BC92", End::Both, b"TTGTGAGTGGAAAGATACAGGACC"),
    ("BC93", End::Both, b"AGTTTCCATCACTTCAGACTTGGG"),
    ("BC94", End::Both, b"GATTGTCCTCAAACTGCCACCTAC"),
    ("BC95", End::Both, b"CCTGTCTGGAAGAAGAATGGACTT"),
    ("BC96", End::Both, b"CTGAACGGTCATAGAGTCCACCAT"),
];
