use std::path::Path;

use anyhow::Result;
use cyto_io::open_file_handle;
use derive_more::AddAssign;
use log::debug;

use crate::ReadMate;

#[derive(Clone, Copy, Debug, serde::Serialize)]
pub struct InputRuntimeStatistics {
    /// The input identity
    pub input_id: usize,
    /// The amount of time taken to map the input (seconds)
    pub elapsed_sec: f64,
}

#[derive(Clone, Copy, Debug, serde::Serialize)]
pub struct LibraryStatistics {
    /// The name of the library type
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
    #[serde(flatten)]
    reads: MappingStatistics,
    libraries: &'a [LibraryStatistics],
    inputs: &'a [InputRuntimeStatistics],
}

pub fn write_statistics<P: AsRef<Path>>(
    outdir: P,
    libstats: &[LibraryStatistics],
    mapstats: MappingStatistics,
    runstats: &[InputRuntimeStatistics],
) -> Result<()> {
    impl_write_statistics(
        outdir.as_ref().join("stats").join("mapping.json"),
        MappingFile {
            reads: mapstats,
            libraries: libstats,
            inputs: runstats,
        },
    )
}

fn impl_write_statistics<P: AsRef<Path>, S: serde::Serialize>(path: P, stat: S) -> Result<()> {
    debug!("Saving statistics to: {}", path.as_ref().display());
    let mut handle = open_file_handle(&path)?;
    serde_json::to_writer_pretty(&mut handle, &stat)?;
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
        let inputs = [InputRuntimeStatistics {
            input_id: 0,
            elapsed_sec: 1.5,
        }];
        let file = MappingFile {
            reads,
            libraries: &[],
            inputs: &inputs,
        };
        let json = serde_json::to_value(file).unwrap();
        assert_eq!(json["total_reads"], 10);
        assert_eq!(json["mapped_reads"], 6);
        assert_eq!(json["unmapped"]["missing_whitelist"], 3);
        assert_eq!(json["libraries"], serde_json::json!([]));
        assert_eq!(json["inputs"][0]["elapsed_sec"], 1.5);
    }
}
