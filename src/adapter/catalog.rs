//! The built-in adapter, primer, barcode and flank catalog.
//!
//! Assembled from primary sources and cross-verified: dorado
//! `adapter_primer_kits.cpp` and `utils/barcode_kits.cpp` (kit-14
//! authoritative), Porechop `porechop/adapters.py` (legacy kits and the 96
//! barcodes), the qcat kit YAMLs, and the NCBI UniVec records for the PacBio
//! sequences. The sequences are vendor-published facts; the catalog itself is
//! whittle's own compilation.
//!
//! Every entry names the kits it belongs to, so a preset selects the subset a
//! library can contain. Every entry is searched on both strands, so a rear
//! sequence stored as the reverse complement of a front sequence is found at
//! either end regardless.
//!
//! Barcode dedup rule: one barcode number is a single shared 24 bp oligo across
//! the native, PCR, and rapid kits, which differ only in flank and orientation
//! (native uses the reverse complement). Only the canonical forward oligo is
//! stored, and reverse-complement search covers the rest. Trimming through the
//! flanks removes the barcode whatever its number.

use super::Role;
use super::preset::Kit;

/// One catalog entry: display name, role, member kits, and sequence.
/// Sequences are uppercase nucleotide codes of at least `MIN_PATTERN_LEN`
/// bases, which `preset::tests::entries_are_valid_nucleotide_sequences` and
/// `preset::tests::entries_meet_the_minimum_pattern_length` enforce.
pub(super) type Entry = (&'static str, Role, &'static [Kit], &'static [u8]);

/// Every kit-14 ligation, rapid and PCR kit carries the ligation Y-adapter.
const KIT14: &[Kit] = &[
    Kit::Lsk114,
    Kit::Rad114,
    Kit::Rbk114,
    Kit::Nbd114,
    Kit::Pcb114,
    Kit::Rpb114,
    Kit::Mab114,
];

/// Kits that attach the rapid adapter by transposase or rapid attachment.
const RAPID: &[Kit] = &[Kit::Rad114, Kit::Rbk114, Kit::Rpb114, Kit::Mab114];

/// Every catalog entry in display order. `preset::build` collapses duplicate
/// sequences.
#[rustfmt::skip]
pub(super) const CATALOG: &[Entry] = &[

    // Ligation Y-adapter, two chemistry generations, both kept.
    // SQK-LSK114 and every kit-14 kit [dorado:adapter_primer_kits.cpp(LSK110)]
    ("LSK114_front", Role::Adapter, KIT14, b"CCTGTACTTCGTTCAGTTACGTATTGC"),
    // [dorado:adapter_primer_kits.cpp(LSK110)]
    ("LSK114_rear", Role::Adapter, KIT14, b"AGCAATACGTAACTGAAC"),
    // SQK-LSK108/109, NSK007 (kit 9/10) [porechop:adapters.py(SQK-NSK007_Y_Top)]
    ("LSK109_front", Role::Adapter, KIT14, b"AATGTACTTCGTTCAGTTACGTATTGCT"),
    // [porechop:adapters.py(SQK-NSK007_Y_Bottom)]
    ("LSK109_rear", Role::Adapter, KIT14, b"GCAATACGTAACTGAACGAAGT"),

    // Rapid adapter, which doubles as the 3' flank of the rapid-barcoding kits.
    // SQK-RAD004/RAD114, ULK114, RBK*, RPB*, MAB114
    // [dorado:adapter_primer_kits.cpp(RAD) / porechop(Rapid_adapter)]
    ("RAD", Role::Adapter, RAPID, b"GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"),

    // Direct RNA.
    // SQK-RNA004, RNA004-XL, DRB004 [dorado:adapter_primer_kits.cpp(RNA004)]
    ("RNA004_rear", Role::Adapter, &[Kit::Rna004], b"GGTTGTTTCTGTTGGTGCTG"),

    // PacBio SMRTbell blunt adapter and C2 sequencing primer
    // [NCBI UniVec NGB00972.1, NGB00973.1].
    ("SMRTbell", Role::Adapter, &[Kit::PacBio], b"ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT"),
    ("PacBio_C2", Role::Primer, &[Kit::PacBio], b"AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA"),

    // Primers: PCR, cDNA, and 10X.
    // SQK-LSK114 cDNA [dorado:adapter_primer_kits.cpp(cDNA)]
    ("cDNA_front", Role::Primer, &[Kit::Lsk114], b"TTTCTGTTGGTGCTGATATTGCTGGG"),
    // [dorado:adapter_primer_kits.cpp(cDNA)]
    ("cDNA_rear", Role::Primer, &[Kit::Lsk114], b"ACTTGCCTGTCGCTCTATCTTCTTT"),
    // SQK-PCS114, PCB114 [dorado:adapter_primer_kits.cpp(PCS110)]
    ("PCS110_front", Role::Primer, &[Kit::Pcb114], b"TTTCTGTTGGTGCTGATATTGCTTT"),
    // [dorado:adapter_primer_kits.cpp(PCS110)]
    ("PCS110_rear", Role::Primer, &[Kit::Pcb114], b"ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAGTTTT"),
    // SQK-LSK114 (10X) [dorado:adapter_primer_kits.cpp(GEN10X)]
    ("GEN10X_front", Role::Primer, &[Kit::Lsk114], b"CTACACGACGCTCTTCCGATCT"),
    // [dorado:adapter_primer_kits.cpp(GEN10X)]
    ("GEN10X_rear", Role::Primer, &[Kit::Lsk114], b"GTACTCTGCGTTGATACCACTGCTT"),
    // PCR/cDNA (legacy) [porechop:adapters.py(PCR_1)]; also LWB_1st_REAR.
    ("PCR1_front", Role::Primer, &[Kit::Pcb114, Kit::Rpb114], b"ACTTGCCTGTCGCTCTATCTTC"),
    // [porechop:adapters.py(PCR_1)]
    ("PCR1_rear", Role::Primer, &[Kit::Pcb114, Kit::Rpb114], b"GAAGATAGAGCGACAGGCAAGT"),
    // PCR/cDNA (legacy) [porechop:adapters.py(PCR_2)]; also LWB_2nd_REAR.
    ("PCR2_front", Role::Primer, &[Kit::Pcb114, Kit::Rpb114], b"TTTCTGTTGGTGCTGATATTGC"),
    // [porechop:adapters.py(PCR_2)]
    ("PCR2_rear", Role::Primer, &[Kit::Pcb114, Kit::Rpb114], b"GCAATATCAGCACCAACAGAAA"),
    // cDNA SSP (legacy) [porechop:adapters.py(cDNA_SSP)]
    ("cDNA_SSP", Role::Primer, &[Kit::Pcb114], b"TTTCTGTTGGTGCTGATATTGCTGCCATTACGGCCGGG"),

    // Universal marker-gene primers of the microbial amplicon kit, IUPAC
    // degenerate: 16S 27F and 1492R, ITS1F and ITS4.
    ("16S_27F", Role::Primer, &[Kit::Mab114], b"AGAGTTTGATYMTGGCTCAG"),
    ("16S_1492R", Role::Primer, &[Kit::Mab114], b"TACGGYTACCTTGTTACGACTT"),
    ("ITS1F", Role::Primer, &[Kit::Mab114], b"CTTGGTCATTTAGAGGAAGTAA"),
    ("ITS4", Role::Primer, &[Kit::Mab114], b"TCCTCCGCTTATTGATATGC"),

    // Barcode flanks (dorado kit-14 constants). Trimming through a flank removes the
    // barcode regardless of its number. Flanks shorter than `MIN_PATTERN_LEN`
    // (dorado's NB_1st_REAR, BC_1st_FRONT, RBK_FRONT and RLB_FRONT, 7 to 8 bp)
    // are not searched standalone. NB_1st_REAR is the only one on the insert
    // side of its barcode and is searched within the native barcode construct.
    // Native barcoding (NBD*) [dorado:barcode_kits.cpp(NB_1st_FRONT)]
    ("NB_front", Role::Barcode, &[Kit::Nbd114], b"ATTGCTAAGGTTAA"),
    // Native barcode construct: NB_1st_FRONT, the 24 bp barcode as N, and the
    // 8 bp NB_1st_REAR, so a hit trims through the inner flank.
    // [dorado:barcode_kits.cpp(NB_1st_FRONT, NB_1st_REAR)]
    ("NB_construct", Role::Barcode, &[Kit::Nbd114], b"ATTGCTAAGGTTAANNNNNNNNNNNNNNNNNNNNNNNNCAGCACCT"),
    // PCR barcoding (PBC/BC*) [dorado:barcode_kits.cpp(BC_1st_REAR)]
    ("PBC_rear", Role::Barcode, &[Kit::Pcb114, Kit::Rpb114], b"TTAACCTTTCTGTTGGTGCTGATATTGC"),
    // Rapid barcoding v4 / kit 14 (RBK*) and MAB114 [dorado:barcode_kits.cpp(RBK4_FRONT, MAB_FRONT)]
    ("RBK4_front", Role::Barcode, &[Kit::Rbk114, Kit::Mab114], b"GCTTGGGTGTTTAACC"),
    // Rapid ligation barcoding (RLB, RPB) [dorado:barcode_kits.cpp(RLB_REAR)]
    ("RLB_rear", Role::Barcode, &[Kit::Rpb114], b"CGTTTTTCGTGCGCCGCTTC"),
    // Microbial amplicon barcoding (MAB114) [dorado:barcode_kits.cpp(MAB_REAR)]
    ("MAB_rear", Role::Barcode, &[Kit::Mab114], b"CCATATCCGTGTCGCCCTT"),

    // The 96 canonical forward barcode oligos, from Porechop `adapters.py` and confirmed
    // present in dorado. BC01 to BC24 serve the 24-plex PCR and rapid PCR kits as well.
    ("BC01", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"AAGAAAGTTGTCGGTGTCTTTGTG"),
    ("BC02", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"TCGATTCCGTTTGTAGTCGTCTGT"),
    ("BC03", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"GAGTCTTGTGTCCCAGTTACCAGG"),
    ("BC04", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"TTCGGATTCTATCGTGTTTCCCTA"),
    ("BC05", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"CTTGTCCAGGGTTTGTGTAACCTT"),
    ("BC06", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"TTCTCGCAAAGGCAGAAAGTAGTC"),
    ("BC07", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"GTGTTACCGTGGGAATGAATCCTT"),
    ("BC08", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"TTCAGGGAACAAACCAAGTTACGT"),
    ("BC09", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"AACTAGGCACAGCGAGTCTTGGTT"),
    ("BC10", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"AAGCGTTGAAACCTTTGTCCTCTC"),
    ("BC11", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"GTTTCATCTATCGGAGGGAATGGA"),
    ("BC12", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"CAGGTAGAAAGAAGCAGAATCGGA"),
    ("BC13", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"AGAACGACTTCCATACTCGTGTGA"),
    ("BC14", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"AACGAGTCTCTTGGGACCCATAGA"),
    ("BC15", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"AGGTCTACCTCGCTAACACCACTG"),
    ("BC16", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"CGTCAACTGACAGTGGTTCGTACT"),
    ("BC17", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"ACCCTCCAGGAAAGTACCTCTGAT"),
    ("BC18", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"CCAAACCCAACAACCTAGATAGGC"),
    ("BC19", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"GTTCCTCGTGCAGTGTCAAGAGAT"),
    ("BC20", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"TTGCGTCCTGTTACGAGAACTCAT"),
    ("BC21", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"GAGCCTCTCATTGTCCGTTCTCTA"),
    ("BC22", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"ACCACTGCCATGTATCAAAGTACG"),
    ("BC23", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"CTTACTACCCAGTGAACCTCCTCG"),
    ("BC24", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114, Kit::Pcb114, Kit::Rpb114], b"GCATAGTTCTGCATGATGGGTTAG"),
    ("BC25", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GTAAGTTGGGTATGCAACGCAATG"),
    ("BC26", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CATACAGCGACTACGCATTCTCAT"),
    ("BC27", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CGACGGTTAGATTCACCTCTTACA"),
    ("BC28", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TGAAACCTAAGAAGGCACCGTATC"),
    ("BC29", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CTAGACACCTTGGGTTGACAGACC"),
    ("BC30", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TCAGTGAGGATCTACTTCGACCCA"),
    ("BC31", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TGCGTACAGCAATCAGTTACATTG"),
    ("BC32", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CCAGTAGAAGTCCGACAACGTCAT"),
    ("BC33", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CAGACTTGGTACGGTTGGGTAACT"),
    ("BC34", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GGACGAAGAACTCAAGTCAAAGGC"),
    ("BC35", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CTACTTACGAAGCTGAGGGACTGC"),
    ("BC36", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"ATGTCCCAGTTAGAGGAGGAAACA"),
    ("BC37", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GCTTGCGATTGATGCTTAGTATCA"),
    ("BC38", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"ACCACAGGAGGACGATACAGAGAA"),
    ("BC39", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CCACAGTGTCAACTAGAGCCTCTC"),
    ("BC40", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TAGTTTGGATGACCAAGGATAGCC"),
    ("BC41", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GGAGTTCGTCCAGAGAAGTACACG"),
    ("BC42", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CTACGTGTAAGGCATACCTGCCAG"),
    ("BC43", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CTTTCGTTGTTGACTCGACGGTAG"),
    ("BC44", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"AGTAGAAAGGGTTCCTTCCCACTC"),
    ("BC45", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GATCCAACAGAGATGCCTTCAGTG"),
    ("BC46", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GCTGTGTTCCACTTCATTCTCCTG"),
    ("BC47", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GTGCAACTTTCCCACAGGTAGTTC"),
    ("BC48", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CATCTGGAACGTGGTACACCTGTA"),
    ("BC49", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"ACTGGTGCAGCTTTGAACATCTAG"),
    ("BC50", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"ATGGACTTTGGTAACTTCCTGCGT"),
    ("BC51", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GTTGAATGAGCCTACTGGGTCCTC"),
    ("BC52", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TGAGAGACAAGATTGTTCGTGGAC"),
    ("BC53", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"AGATTCAGACCGTCTCATGCAAAG"),
    ("BC54", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CAAGAGCTTTGACTAAGGAGCATG"),
    ("BC55", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TGGAAGATGAGACCCTGATCTACG"),
    ("BC56", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TCACTACTCAACAGGTGGCATGAA"),
    ("BC57", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GCTAGGTCAATCTCCTTCGGAAGT"),
    ("BC58", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CAGGTTACTCCTCCGTGAGTCTGA"),
    ("BC59", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TCAATCAAGAAGGGAAAGCAAGGT"),
    ("BC60", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CATGTTCAACCAAGGCTTCTATGG"),
    ("BC61", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"AGAGGGTACTATGTGCCTCAGCAC"),
    ("BC62", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CACCCACACTTACTTCAGGACGTA"),
    ("BC63", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TTCTGAAGTTCCTGGGTCTTGAAC"),
    ("BC64", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GACAGACACCGTTCATCGACTTTC"),
    ("BC65", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TTCTCAGTCTTCCTCCAGACAAGG"),
    ("BC66", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CCGATCCTTGTGGCTTCTAACTTC"),
    ("BC67", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GTTTGTCATACTCGTGTGCTCACC"),
    ("BC68", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GAATCTAAGCAAACACGAAGGTGG"),
    ("BC69", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TACAGTCCGAGCCTCATGTGATCT"),
    ("BC70", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"ACCGAGATCCTACGAATGGAGTGT"),
    ("BC71", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CCTGGGAGCATCAGGTAGTAACAG"),
    ("BC72", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TAGCTGACTGTCTTCCATACCGAC"),
    ("BC73", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"AAGAAACAGGATGACAGAACCCTC"),
    ("BC74", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TACAAGCATCCCAACACTTCCACT"),
    ("BC75", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GACCATTGTGATGAACCCTGTTGT"),
    ("BC76", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"ATGCTTGTTACATCAACCCTGGAC"),
    ("BC77", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CGACCTGTTTCTCAGGGATACAAC"),
    ("BC78", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"AACAACCGAACCTTTGAATCAGAA"),
    ("BC79", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TCTCGGAGATAGTTCTCACTGCTG"),
    ("BC80", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CGGATGAACATAGGATAGCGATTC"),
    ("BC81", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CCTCATCTTGTGAAGTTGTTTCGG"),
    ("BC82", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"ACGGTATGTCGAGTTCCAGGACTA"),
    ("BC83", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TGGCTTGATCTAGGTAAGGTCGAA"),
    ("BC84", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GTAGTGGACCTAGAACCTGTGCCA"),
    ("BC85", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"AACGGAGGAGTTAGTTGGATGATC"),
    ("BC86", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"AGGTGATCCCAACAAGCGTAAGTA"),
    ("BC87", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TACATGCTCCTGTTGTTAGGGAGG"),
    ("BC88", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TCTTCTACTACCGATCCGAAGCAG"),
    ("BC89", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"ACAGCATCAATGTTTGGCTAGTTG"),
    ("BC90", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GATGTAGAGGGTACGGTTTGAGGC"),
    ("BC91", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GGCTCCATAGGAACTCACGCTACT"),
    ("BC92", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"TTGTGAGTGGAAAGATACAGGACC"),
    ("BC93", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"AGTTTCCATCACTTCAGACTTGGG"),
    ("BC94", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"GATTGTCCTCAAACTGCCACCTAC"),
    ("BC95", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CCTGTCTGGAAGAAGAATGGACTT"),
    ("BC96", Role::Barcode, &[Kit::Rbk114, Kit::Nbd114], b"CTGAACGGTCATAGAGTCCACCAT"),

    // The 24 barcodes of the microbial amplicon kit [dorado:barcode_kits.cpp(TP01..TP24)].
    ("TP01", Role::Barcode, &[Kit::Mab114], b"GCACCTGGAACTTGTGCCTTCCAC"),
    ("TP02", Role::Barcode, &[Kit::Mab114], b"CCGAAATAGGTTATCTGTTGTTGT"),
    ("TP03", Role::Barcode, &[Kit::Mab114], b"ATCAATCGCTGGACGATGGATTAG"),
    ("TP04", Role::Barcode, &[Kit::Mab114], b"CCACCCGCTCCTGCCGGTGGGCGT"),
    ("TP05", Role::Barcode, &[Kit::Mab114], b"AGACTCTTGGGCTCGCCACGTCCC"),
    ("TP06", Role::Barcode, &[Kit::Mab114], b"TCTGTATCCGGAGACGGGATGGAC"),
    ("TP07", Role::Barcode, &[Kit::Mab114], b"TTTCGGATCAATCGACCGCAAACG"),
    ("TP08", Role::Barcode, &[Kit::Mab114], b"ACTCAAACATTCTGTTAGATCGCG"),
    ("TP09", Role::Barcode, &[Kit::Mab114], b"AAATGGAACCCGGATATGTTTACT"),
    ("TP10", Role::Barcode, &[Kit::Mab114], b"TAAATCGACCTATGATGAACACAG"),
    ("TP11", Role::Barcode, &[Kit::Mab114], b"ACATGTTGGAGTGAAAGTCGGGTA"),
    ("TP12", Role::Barcode, &[Kit::Mab114], b"CCTGGACCACGATCATTGTAACAT"),
    ("TP13", Role::Barcode, &[Kit::Mab114], b"TATGGTGGATCTCCCTCTATCTTC"),
    ("TP14", Role::Barcode, &[Kit::Mab114], b"AAGTAAATGGGACGCCCACTCCGA"),
    ("TP15", Role::Barcode, &[Kit::Mab114], b"TGTTCGCGGCTTGATCTAATATTA"),
    ("TP16", Role::Barcode, &[Kit::Mab114], b"AGAGAGCTTCCCGGGAGGGTGGTC"),
    ("TP17", Role::Barcode, &[Kit::Mab114], b"TTGTGAATATCTGTCACAAACACC"),
    ("TP18", Role::Barcode, &[Kit::Mab114], b"CAATCGTACCAGGGAACATAAAGT"),
    ("TP19", Role::Barcode, &[Kit::Mab114], b"CACACCCAAACAATATGGACCCGT"),
    ("TP20", Role::Barcode, &[Kit::Mab114], b"AATAACCACATCCGCCCTCCGCAC"),
    ("TP21", Role::Barcode, &[Kit::Mab114], b"TCCTAATAATGTGTAGATCGGTCC"),
    ("TP22", Role::Barcode, &[Kit::Mab114], b"AGTCGATGGAACAAGAGAAGTTAT"),
    ("TP23", Role::Barcode, &[Kit::Mab114], b"AAACTCACTGTATGTCGTTTCTAT"),
    ("TP24", Role::Barcode, &[Kit::Mab114], b"TGACATCACTGATCGAGGAAGATC"),
];
