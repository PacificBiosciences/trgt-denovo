use crate::allele::Allele;

pub fn get_per_allele_reads(alleles: &[Allele]) -> Vec<usize> {
    alleles.iter().map(|a| a.read_aligns.len()).collect()
}

pub fn get_total_reads(alleles: &[Allele]) -> usize {
    get_per_allele_reads(alleles).iter().sum()
}

pub fn get_dropout_prob(alleles: &[Allele]) -> f64 {
    assert!(!alleles.is_empty());
    let allele_frac: f64 = 0.5;
    let num_reads = get_total_reads(alleles) as f64;
    allele_frac.powf(num_reads)
}

pub fn get_allele_freqs(alleles: &[Allele]) -> Vec<f64> {
    let num_reads = get_total_reads(alleles) as f64;
    let allele_freqs: Vec<f64> = alleles
        .iter()
        .map(|a| a.read_aligns.len() as f64 / num_reads)
        .collect();
    allele_freqs
}

#[cfg(test)]
mod tests {}
