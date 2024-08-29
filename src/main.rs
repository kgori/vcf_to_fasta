#![allow(unused_variables)]
use clap::Parser;

use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use rust_htslib::faidx;

use itertools::Itertools;

use std::collections::HashMap;

use bio::io::fasta::{Record, Writer};
use std::fs;
use std::io;
use std::io::Write;

#[derive(Parser)]
struct Cli {
    #[clap(parse(from_os_str))]
    fasta_file: std::path::PathBuf,

    #[clap(parse(from_os_str))]
    vcf_file: std::path::PathBuf,

    #[clap(parse(from_os_str))]
    out_file: std::path::PathBuf,

    #[clap(short = 'a', long, default_value_t = 3)]
    min_alt_reads: i32,

    #[clap(short = 'm', long, default_value_t = 0.1)]
    min_vaf: f64,

    #[clap(short = 'f', long)]
    write_fasta: bool,

    #[clap(short = 'l', long)]
    write_lists: bool,
}

#[derive(Debug)]
#[allow(dead_code)]
struct Variant {
    ref_allele: String,
    alt_allele: String,
    chrom: String,
    pos: i64,
    vaf: f64,
}

impl Variant {
    fn new(ref_a: &String, alt_a: &String, chr: &String, position: i64, vaf: f64) -> Self {
        Self {
            ref_allele: ref_a.to_owned(),
            alt_allele: alt_a.to_owned(),
            chrom: chr.to_owned(),
            pos: position,
            vaf: vaf,
        }
    }
}

fn main() -> Result<(), anyhow::Error> {
    let args = Cli::parse();

    let fasta = faidx::Reader::from_path(args.fasta_file)?;
    let mut vcf = bcf::Reader::from_path(args.vcf_file)?;
    let header = vcf.header().clone();
    let sample_count: usize = header.sample_count().try_into()?;
    let sample_names = header
        .samples()
        .iter()
        .map(|name| std::string::String::from(std::str::from_utf8(name).unwrap()))
        .collect_vec();

    let mut variant_map = HashMap::<String, Vec<Variant>>::new();
    for name in &sample_names {
        variant_map.insert(name.to_string(), Vec::new());
    }

    let seq = fasta.fetch_seq_string("chrM", 0, 16726)?;
    let seq = seq.to_uppercase();
    let seq = seq.as_bytes();

    for record in vcf.records() {
        let record = record?;
        let rid = record.rid().unwrap();
        let chrom = header.rid2name(rid)?;
        let chrom = std::string::String::from(std::str::from_utf8(chrom).unwrap());
        let pos = record.pos();
        let (ref_allele, alt_allele) = record
            .alleles()
            .iter()
            .map(|a| std::string::String::from(std::str::from_utf8(a).unwrap()))
            .collect_tuple()
            .unwrap();

        let read_depths = record
            .format(b"NR")
            .integer()
            .expect("Couldn't retrieve NR field");
        let variant_depths = record
            .format(b"NV")
            .integer()
            .expect("Couldn't retrieve NV field");

        // println!("Chrom:{} Pos:{} Ref:{} Alt:{}", chrom, pos, ref_allele, alt_allele);
        for sample_id in 0..sample_count {
            let read_depth = *read_depths[sample_id].get(0).unwrap() as f64;
            let variant_depth = *variant_depths[sample_id].get(0).unwrap() as f64;
            let vaf = variant_depth / read_depth;
            let sample_name = &sample_names[sample_id];

            if vaf > args.min_vaf && variant_depth >= args.min_alt_reads.into() {
                let new_variant = Variant::new(&ref_allele, &alt_allele, &chrom, pos, vaf);
                variant_map
                    .entry(sample_name.to_string())
                    .or_insert(Vec::new())
                    .push(new_variant);
            }
        }
    }
    println!("variant_map size = {}", variant_map.len());

    if args.write_fasta {
        let path = args.out_file.clone();
        let file = fs::File::create(path).unwrap();
        let handle = io::BufWriter::new(file);
        let mut writer = Writer::new(handle);

        for (name, variants) in variant_map.iter() {
            let mutable_seq: &mut [u8] = &mut seq.to_owned();
            for variant in variants {
                // Skip hyper-variable region
                if variant.pos < 16129 || variant.pos > 16430 {
                    let seq_at_pos: u8 = seq[variant.pos as usize];
                    let ref_char: u8 = *variant.ref_allele.as_bytes().get(0).unwrap();
                    let alt_char: u8 = *variant.alt_allele.as_bytes().get(0).unwrap();
                    if seq_at_pos != ref_char {
                        println!(
                            "Pos:{}, Seq:{:?} != VcfRef:{:?} (VcfAlt:{:?})",
                            variant.pos,
                            std::string::String::from_utf8((&[seq_at_pos]).to_vec()),
                            std::string::String::from_utf8((&[ref_char]).to_vec()),
                            std::string::String::from_utf8((&[alt_char]).to_vec()),
                        );
                    }
                    assert!(seq_at_pos == ref_char);
                    mutable_seq[variant.pos as usize] = alt_char;
                }
            }
            let fasta_record = Record::with_attrs(name, None, mutable_seq);
            let write_result = writer.write_record(&fasta_record);
            assert!(write_result.is_ok());
        }
    }

    if args.write_lists {
        for (name, variants) in variant_map {
            let mut path = args.out_file.clone();
            path.set_extension("listsdir");
            std::fs::create_dir_all(&path)?;
            let outfile = path.join(name + ".txt");
            println!("outfile = {:?}", outfile);
            let file = fs::File::create(outfile).unwrap();
            let mut writer = io::BufWriter::new(file);
            for variant in variants {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}",
                    variant.chrom, variant.pos, variant.ref_allele, variant.alt_allele, variant.vaf
                )?;
            }
        }
    }

    Ok(())
}
