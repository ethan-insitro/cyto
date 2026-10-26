use std::path::Path;

use anyhow::Result;
use cyto_io::open_file_handle;
use derive_more::AddAssign;
use log::debug;

use crate::ReadMate;

#[derive(Clone, Debug, serde::Serialize)]
pub struct InputRuntimeStatistics {
    /// The input files, as given on the command line
    pub inputs: Vec<String>,
    /// The amount of time taken to map all inputs (seconds)
    pub elapsed_sec: f64,
    /// The time taken to map each input, in `inputs` order. `None` for FASTQ, whose
    /// R1/R2 pairs are all mapped in one pass
    pub per_input_elapsed_sec: Option<Vec<f64>>,
}

#[derive(Clone, Copy, Debug, serde::Serialize)]
pub struct LibraryStatistics {
    /// The name of the library type: its key in `mapping.json`'s `libraries`
    #[serde(skip_serializing)]
    pub name: &'static str,

    /// The total number of elements in the library
    pub total_elem: usize,

    /// The number of aggregated elements in the library
    ///
    /// This can vary in meaning depending on the library type.
    /// For GEX - it's the number of genes
    /// For Probes - it's the number of demultiplexing probes (not their sequences)
    pub total_aggr: usize,

    /// The total number of hash keys in this library
    pub total_hash: usize,

    /// The mapped position of this library in the query constructs
    pub position: usize,

    /// The read mate searched for this library
    pub mate: ReadMate,

    /// The remap window size used for this library
    pub window: usize,

    /// Whether this library used exact matching
    pub exact: bool,

    /// The amount of time taken to initialize the library (seconds)
    pub init_time: f64,
}

/// Aggregate mapping statistics for a library
#[derive(Clone, Copy, Default, Debug, AddAssign, serde::Serialize)]
pub struct MappingStatistics {
    pub total_reads: usize,
    pub mapped_reads: usize,
    pub unmapped: UnmappedStatistics,
}
impl MappingStatistics {
    /// Calculate the percentage mapped
    pub fn frac_mapped(&self) -> f64 {
        self.mapped_reads as f64 / self.total_reads as f64
    }
    /// Resets all counters on the statistics
    pub fn clear(&mut self) {
        *self = Self::default();
    }
}

/// Reads that failed each mapping check. A read can fail more than one.
#[derive(Clone, Copy, Default, Debug, AddAssign, serde::Serialize)]
pub struct UnmappedStatistics {
    pub missing_probe: usize,
    pub missing_feature: usize,
    pub missing_whitelist: usize,
    pub failed_umi_qual: usize,
    pub umi_truncated: usize,
}

/// `stats/mapping.json`: everything the mapping step measured
#[derive(serde::Serialize)]
struct MappingFile<'a> {
    total_reads: usize,
    mapped_reads: usize,
    fraction_mapped: f64,
    unmapped: UnmappedStatistics,
    #[serde(serialize_with = "by_name")]
    libraries: &'a [LibraryStatistics],
    #[serde(flatten)]
    runtime: &'a InputRuntimeStatistics,
}

impl<'a> MappingFile<'a> {
    fn new(
        reads: MappingStatistics,
        libraries: &'a [LibraryStatistics],
        runtime: &'a InputRuntimeStatistics,
    ) -> Self {
        let fraction_mapped = if reads.total_reads > 0 {
            reads.mapped_reads as f64 / reads.total_reads as f64
        } else {
            0.0
        };
        Self {
            total_reads: reads.total_reads,
            mapped_reads: reads.mapped_reads,
            fraction_mapped,
            unmapped: reads.unmapped,
            libraries,
            runtime,
        }
    }
}

/// Writes the libraries as one object keyed by name, in mapping order
fn by_name<S: serde::Serializer>(
    libs: &&[LibraryStatistics],
    serializer: S,
) -> Result<S::Ok, S::Error> {
    serializer.collect_map(libs.iter().map(|lib| (lib.name, lib)))
}

pub fn write_statistics<P: AsRef<Path>>(
    outdir: P,
    libstats: &[LibraryStatistics],
    mapstats: MappingStatistics,
    runstats: &InputRuntimeStatistics,
) -> Result<()> {
    let path = outdir.as_ref().join("stats").join("mapping.json");
    debug!("Saving statistics to: {}", path.display());
    let mut handle = open_file_handle(&path)?;
    let file = MappingFile::new(mapstats, libstats, runstats);
    serde_json::to_writer_pretty(&mut handle, &file)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mapping_json_holds_the_counts_libraries_and_inputs() {
        let reads = MappingStatistics {
            total_reads: 10,
            mapped_reads: 6,
            unmapped: UnmappedStatistics {
                missing_whitelist: 3,
                ..Default::default()
            },
        };
        let runtime = InputRuntimeStatistics {
            inputs: vec!["a.cbq".to_string(), "b.cbq".to_string()],
            elapsed_sec: 1.5,
            per_input_elapsed_sec: Some(vec![1.0, 0.5]),
        };
        let libraries = [LibraryStatistics {
            name: "gex",
            total_elem: 5,
            total_aggr: 2,
            total_hash: 5,
            position: 0,
            mate: ReadMate::R2,
            window: 1,
            exact: false,
            init_time: 0.5,
        }];
        let file = MappingFile::new(reads, &libraries, &runtime);
        let json = serde_json::to_value(file).unwrap();
        assert_eq!(json["total_reads"], 10);
        assert_eq!(json["mapped_reads"], 6);
        assert_eq!(json["fraction_mapped"], 0.6);
        assert_eq!(json["unmapped"]["missing_whitelist"], 3);
        assert_eq!(json["libraries"]["gex"]["total_aggr"], 2);
        assert!(json["libraries"]["gex"].get("name").is_none());
        assert_eq!(json["inputs"], serde_json::json!(["a.cbq", "b.cbq"]));
        assert_eq!(json["elapsed_sec"], 1.5);
        assert_eq!(json["per_input_elapsed_sec"], serde_json::json!([1.0, 0.5]));

        let fastq = InputRuntimeStatistics {
            per_input_elapsed_sec: None,
            ..runtime
        };
        let empty =
            serde_json::to_value(MappingFile::new(MappingStatistics::default(), &[], &fastq))
                .unwrap();
        assert_eq!(empty["fraction_mapped"], 0.0);
        assert!(empty["per_input_elapsed_sec"].is_null());
    }
}
