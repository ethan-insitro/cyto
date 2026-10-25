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
///
/// Serialized as `MappingMapFile`, the `mapping_map.json` form.
#[derive(Clone, Copy, Default, Debug, AddAssign, serde::Serialize)]
#[serde(into = "MappingMapFile")]
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

/// `stats/mapping_map.json`: fractions are 0 without reads
#[derive(serde::Serialize)]
struct MappingMapFile {
    total_reads: usize,
    mapped_reads: usize,
    unmapped_reads: usize,
    mapped_reads_frac: f64,
    unmapped_reads_frac: f64,
    /// Largest first
    unmapped: Vec<UnmappedReason>,
}

/// One reason reads failed to map
#[derive(serde::Serialize)]
struct UnmappedReason {
    reason: &'static str,
    reads: usize,
    frac_of_reads: f64,
    frac_of_unmapped: f64,
}

impl From<MappingStatistics> for MappingMapFile {
    fn from(m: MappingStatistics) -> Self {
        let unmapped_reads = m.total_reads - m.mapped_reads;
        let frac = |x: usize, total: usize| {
            if total > 0 {
                x as f64 / total as f64
            } else {
                0.0
            }
        };
        let u = m.unmapped;
        let mut reasons = [
            ("missing_probe", u.missing_probe),
            ("missing_feature", u.missing_feature),
            ("missing_whitelist", u.missing_whitelist),
            ("failed_umi_qual", u.failed_umi_qual),
            ("umi_truncated", u.umi_truncated),
        ];
        reasons.sort_by_key(|&(_, reads)| std::cmp::Reverse(reads)); // stable
        Self {
            total_reads: m.total_reads,
            mapped_reads: m.mapped_reads,
            unmapped_reads,
            mapped_reads_frac: frac(m.mapped_reads, m.total_reads),
            unmapped_reads_frac: frac(unmapped_reads, m.total_reads),
            unmapped: reasons
                .into_iter()
                .map(|(reason, reads)| UnmappedReason {
                    reason,
                    reads,
                    frac_of_reads: frac(reads, m.total_reads),
                    frac_of_unmapped: frac(reads, unmapped_reads),
                })
                .collect(),
        }
    }
}

/// Aggregate mapping statistics for a library focused on unmapped reads
#[derive(Clone, Copy, Default, Debug, AddAssign)]
pub struct UnmappedStatistics {
    pub missing_probe: usize,
    pub missing_feature: usize,
    pub missing_whitelist: usize,
    pub failed_umi_qual: usize,
    pub umi_truncated: usize,
}

/// `stats/mapping_lib.json`: the libraries keyed by name, in mapping order
#[derive(serde::Serialize)]
struct MappingLibFile<'a> {
    #[serde(serialize_with = "by_name")]
    libraries: &'a [LibraryStatistics],
}
fn by_name<S: serde::Serializer>(libs: &&[LibraryStatistics], s: S) -> Result<S::Ok, S::Error> {
    s.collect_map(libs.iter().map(|lib| (lib.name, lib)))
}

/// `stats/mapping_run.json`
#[derive(serde::Serialize)]
struct MappingRunFile<'a> {
    inputs: &'a [InputRuntimeStatistics],
}

pub fn write_statistics<P: AsRef<Path>>(
    outdir: P,
    libstats: &[LibraryStatistics],
    mapstats: MappingStatistics,
    runstats: &[InputRuntimeStatistics],
) -> Result<()> {
    let stats_outdir = outdir.as_ref().join("stats");

    impl_write_statistics(
        stats_outdir.join("mapping_lib.json"),
        MappingLibFile {
            libraries: libstats,
        },
    )?;
    impl_write_statistics(stats_outdir.join("mapping_map.json"), mapstats)?;
    impl_write_statistics(
        stats_outdir.join("mapping_run.json"),
        MappingRunFile { inputs: runstats },
    )?;
    Ok(())
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
    fn mapping_map_lists_unmapped_reasons_largest_first() {
        let stats = MappingStatistics {
            total_reads: 10,
            mapped_reads: 6,
            unmapped: UnmappedStatistics {
                missing_feature: 1,
                missing_whitelist: 3,
                ..Default::default()
            },
        };
        let json = serde_json::to_value(stats).unwrap();
        assert_eq!(json["unmapped_reads"], 4);
        assert_eq!(json["unmapped"][0]["reason"], "missing_whitelist");
        assert_eq!(json["unmapped"][0]["frac_of_reads"], 0.3);
        assert_eq!(json["unmapped"][0]["frac_of_unmapped"], 0.75);
        assert_eq!(json["unmapped"][1]["reason"], "missing_feature");
        assert_eq!(json["unmapped"].as_array().unwrap().len(), 5);
        assert_eq!(
            serde_json::to_value(MappingStatistics::default()).unwrap()["mapped_reads_frac"],
            0.0
        );
    }
}
