var previousAlertMultipleSeqs = false;
var previousAlertCharacters = false;

const codonDict = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
};

// Calculate codon usage from a DNA sequence
function calculateCodonUsage(sequence) {
    var codonCount = {};
    var aminoAcidCount = {};
    var totalCodons = 0;

    for (const codon in codonDict) {
        codonCount[codon] = 0;
        aminoAcidCount[codonDict[codon]] = 0;
    }

    for (let i = 0; i < sequence.length; i += 3) {
        const codon = sequence.slice(i, i + 3);
        if (codon in codonDict) {
            codonCount[codon] += 1;
            aminoAcidCount[codonDict[codon]] += 1;
            totalCodons += 1
        }
    }

    var results_table = document.getElementById("results-table")
    results_table.innerHTML = "";

    var thead = results_table.createTHead();
    thead.classList.add("table-dark");
    var row = thead.insertRow(0);

    for (const val of ["Codon", "Amino Acid", "Count", "Total Usage", "Relative Usage"]) {
        cell = row.insertCell();
        cell.innerHTML = val;
    }

    var tbody = results_table.createTBody();

    for (const codon in codonCount) {
        amino_acid = codonDict[codon];
        result_count = codonCount[codon];
        result_total_percent = ((result_count / totalCodons) * 100);
        result_relative_percent = ((result_count / aminoAcidCount[amino_acid]) * 100);

        if (isNaN(result_total_percent)) {
            result_total_percent = 0;
        }
        
        if (isNaN(result_relative_percent)) {
            result_relative_percent = 0;
        }

        var row = tbody.insertRow();

        for (const val of [codon, amino_acid, result_count, result_total_percent.toFixed(3), result_relative_percent.toFixed(3)]) {
            var cell = row.insertCell();
            cell.innerHTML = val;
        }   
    }

    $("#results-card").show(500)

    return codonCount;
}

$("#button-calculate").click(function() {
    var sequence = "";
    var n_seqs = 0;

    var lines = $("#sequence").val().toUpperCase().replaceAll("U", "T").split("\n");
    for (var i = 0; i < lines.length; i++) {
        if (lines[i].trim()[0] != ">") {
            sequence += lines[i].trim();
        }
        else {
            n_seqs += 1;
        }
    }

    if (n_seqs > 1 && previousAlertMultipleSeqs == false) {
        $("#alert-placeholder").append('<div id="alert" class="alert alert-warning alert-dismissible fade show" role="alert"><strong>Note:</strong> You have entered multiple sequences. This tool concatenates all input sequences together. You should ensure all of your sequences are in the correct reading frame, etc...<button type="button" class="btn-close" data-bs-dismiss="alert" aria-label="Close"></button></div>');
        previousAlertMultipleSeqs = true;
    }

    if (/[^ACGT]/.test(sequence) && previousAlertCharacters == false) {
        $("#alert-placeholder").append('<div id="alert" class="alert alert-warning alert-dismissible fade show" role="alert"><strong>Note:</strong> Your input text contains characters other than A, C, G, T, and U. Check your input.<button type="button" class="btn-close" data-bs-dismiss="alert" aria-label="Close"></button></div>');
        previousAlertCharacters = true;
    }

    calculateCodonUsage(sequence);
})

$("#button-clear").click(function() {
    $("#sequence").val("");
    $("#results-card").hide(500)
    $("#results").text("");
})