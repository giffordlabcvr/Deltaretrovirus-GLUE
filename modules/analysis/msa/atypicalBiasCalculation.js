/*
 * ====================================================================================
 * Codon Bias Conservation Analysis
 * ------------------------------------------------------------------------------------
 * GLUE Module Type: ecmaFunctionInvoker
 *
 * This function analyzes codon usage in protein-coding features to identify 
 * evolutionarily conserved codons that deviate from the dominant compositional 
 * bias (e.g. A-rich, C-rich) of the viral genome. These codons may indicate 
 * underlying functional constraints at the nucleotide level, such as splicing 
 * elements or RNA structures.
 *
 * ------------------------------------------------------------------------------------
 * FUNCTION SIGNATURE:
 *
 *     calculateBias(alignmentName, refseqName, featureName,
 *                   biasType, aaConservationCutoff,
 *                   codonConservationCutoff, outputMode)
 *
 * ------------------------------------------------------------------------------------
 * PARAMETERS:
 *
 * alignmentName            [String]   Name of the alignment (e.g. 'AL_TREE_EIAV_Eu1')
 *
 * refseqName               [String]   Name of the reference sequence (e.g. 'REF_EIAV_Eu1_JX480632')
 *
 * featureName              [String]   Coding feature to analyze (e.g. 'gag', 'pol')
 *
 * biasType                 [String]   Expected compositional bias for the virus group.
 *                                     Determines which codons are considered 'atypical'.
 *                                     Accepted values:
 *                                       - "A-rich"   (default for lentiviruses)
 *                                       - "C-rich"   (e.g. deltaretroviruses)
 *                                       - "G-rich"
 *                                       - "T-rich"
 *
 * aaConservationCutoff     [Number]   Minimum % of sequences with the same amino acid 
 *                                     at a codon position to consider it 'AA conserved'.
 *                                     Default: 90
 *
 * codonConservationCutoff  [Number]   Minimum % of sequences with the same codon among 
 *                                     AA-conserved sequences to evaluate atypicality.
 *                                     Default: 98
 *
 * outputMode               [String]   Controls which results are returned:
 *                                       - "all"       (default) – all codon positions
 *                                       - "atypical"  – only conserved, atypical codons
 *
 * ------------------------------------------------------------------------------------
 * RETURNS:
 *
 *   A table with the following columns:
 *     - Gene              → feature name
 *     - Codon Label       → position within feature (1-based)
 *     - Genome Nt Pos     → genome-level coordinate of first nucleotide in codon
 *     - Codon             → majority codon at position (if defined)
 *     - Codon %           → % of sequences with the majority codon
 *     - Amino Acid        → consensus amino acid
 *     - AA %              → % of sequences encoding the consensus AA
 *     - Atypical?         → whether the majority codon is atypical
 *     - Depth             → number of valid codons at this position
 *     - Outcome           → interpretive description of conservation pattern
 *
 * ------------------------------------------------------------------------------------
 * EXAMPLE INVOCATION (GLUE COMMAND LINE):
 *
 *     GLUE> module atypicalBiasCalculation invoke-function calculateBias \
 *             AL_TREE_EIAV_Eu1 REF_EIAV_Eu1_JX480632 gag A-rich 90 98 all
 *
 * ------------------------------------------------------------------------------------
 * NOTES:
 *
 * - For six-fold degenerate amino acids (L, S, R), conservation is assessed
 *   within codon subsets to avoid confounding from functional switching between
 *   disjoint synonymous sets.
 *
 * - Feature coordinates must be mapped on the reference sequence using GLUE's
 *   feature-location system.
 *
 * - The script requires segments have been generated in the specified alignment, and
 *   sufficient depth to support analysis.
 *
 * ====================================================================================
 */


function calculateBias(alignmentName, refseqName, featureName, biasType, aaConservationCutoff, codonConservationCutoff, outputMode) {

    biasType = biasType || "A-rich";
    aaConservationCutoff = aaConservationCutoff || 90;
    codonConservationCutoff = codonConservationCutoff || 98;
    outputMode = outputMode || "all";  // default: output all positions

    var results = [];
    glue.log("INFO", "Starting codon bias analysis for alignment: " + alignmentName);

    // Fetch codon metadata from reference
    var codonMetadataMap = getCodonMetadata(refseqName, featureName);

    // Step 1: Fetch amino acid frequencies from the alignment
    var aaFrequencyResult;
    glue.inMode("alignment/" + alignmentName, function() {
        aaFrequencyResult = glue.tableToObjects(
            glue.command(["amino-acid", "frequency", "-r", refseqName, "-f", featureName])
        );
    });

    // Step 2: Process each codon position
    _.each(aaFrequencyResult, function(row) {
        var codonPosition = row.codon;
        var consensusAA = row.aminoAcid;
        var pctMembers = parseFloat(row.pctMembers);

        var codonMeta = codonMetadataMap[codonPosition];

        var result = {
            Gene: featureName,
            "Codon Label": codonPosition,
            "Genome Nt Pos": codonMeta ? codonMeta.refNt : "N/A",
            "Codon": "",
            "Codon %": "",
            "Amino Acid": consensusAA,
            "AA %": pctMembers.toFixed(2),
            "Atypical?": "",
            "Depth": "",
            "Outcome": ""
        };

        if (pctMembers >= aaConservationCutoff) {
            var expected = getExpectedCodons(consensusAA);
            if (!expected || expected.degeneracy === 1) {
                result["Outcome"] = "Skipped: single-codon AA";
            } else {
                var codonCounts = countCodonsForAA(alignmentName, refseqName, featureName, codonPosition, expected.codons);
                var total = _.reduce(codonCounts, function(sum, count) { return sum + count; }, 0);

                result["Depth"] = total;

                if (total === 0) {
                    result["Codon"] = "N/A";
                    result["Codon %"] = "0.0";
                    result["Atypical?"] = "N/A";
                    result["Outcome"] = "No valid codons found";
                } else if (expected.sets && expected.degeneracy === 6) {
                    // Special handling for six-fold degenerate codons
                    var setFrequencies = expected.sets.map(function(set) {
                        var count = 0;
                        _.each(set, function(codon) {
                            count += codonCounts[codon] || 0;
                        });
                        return count;
                    });

                    var dominantSetIndex = _.indexOf(setFrequencies, _.max(setFrequencies));
                    var dominantSet = expected.sets[dominantSetIndex];

                    var dominantCodon = null;
                    var dominantCodonFreq = 0;

                    _.each(dominantSet, function(codon) {
                        var count = codonCounts[codon] || 0;
                        var freq = (count / total) * 100;
                        if (freq > dominantCodonFreq) {
                            dominantCodon = codon;
                            dominantCodonFreq = freq;
                        }
                    });

                    result["Codon"] = dominantCodon;
                    result["Codon %"] = dominantCodonFreq.toFixed(2);
                    var isAtyp = isAtypical(dominantCodon, consensusAA, biasType);
                    result["Atypical?"] = isAtyp ? "Yes" : "No";

                    if (dominantCodonFreq >= codonConservationCutoff) {
                        result["Outcome"] = "Six-fold codon conserved (" + (isAtyp ? "atypical" : "typical") + ")";
                    } else {
                        result["Outcome"] = "Six-fold codon conservation < cutoff";
                    }
                } else {
                    // Standard codon handling
                    var dominantCodon = null;
                    var dominantCodonFreq = 0;

                    _.each(codonCounts, function(count, codon) {
                        var freq = (count / total) * 100;
                        if (freq > dominantCodonFreq) {
                            dominantCodon = codon;
                            dominantCodonFreq = freq;
                        }
                    });

                    result["Codon"] = dominantCodon;
                    result["Codon %"] = dominantCodonFreq.toFixed(2);
                    var isAtyp = isAtypical(dominantCodon, consensusAA, biasType);
                    result["Atypical?"] = isAtyp ? "Yes" : "No";

                    if (dominantCodonFreq >= codonConservationCutoff) {
                        result["Outcome"] = "Single " + (isAtyp ? "atypical" : "typical") + " codon conserved";
                    } else {
                        result["Outcome"] = "Codon conservation < codonConservationCutoff%";
                    }
                }
            }
        } else {
            result["Outcome"] = "Consensus AA: " + consensusAA + " (polymorphic)";
        }

        // Output filtering
        var include = true;
        if (outputMode === "atypical") {
            include = result["Outcome"] && result["Outcome"].indexOf("atypical") !== -1;
        }

        if (include) {
            results.push(result);
        }
    });

    glue.log("INFO", "Codon bias analysis complete.");
    return results;
}


// Utility: Build codon metadata map from reference feature
function getCodonMetadata(refseqName, featureName) {
    var codonMap = {};

    glue.inMode("reference/" + refseqName + "/feature-location/" + featureName, function() {
        var result = glue.tableToObjects(
            glue.command(["amino-acid"])
        );

        _.each(result, function(row) {
            codonMap[row.codonLabel] = {
                refNt: row.refNt,
                codonNts: row.codonNts,
                aminoAcid: row.aminoAcid
            };
        });
    });

    return codonMap;
}

// Utility: Export alignment column at a codon position
function exportAlignmentColumn(alignmentName, refseqName, featureName, position) {
    var exportResult;
    glue.inMode("module/fastaAlignmentExporter", function() {
        exportResult = glue.command([
            "export", alignmentName,
            "-r", refseqName,
            "-f", featureName,
            "-l", position, position,
            "-a", "-p"
        ]);
    });
    return exportResult.nucleotideFasta.sequences;
}

// Utility: Count codons only for those encoding the conserved amino acid
function countCodonsForAA(alignmentName, refseqName, featureName, position, expectedCodons) {
    var sequences = exportAlignmentColumn(alignmentName, refseqName, featureName, position);
    var codonCounts = {};
    _.each(expectedCodons, function(codon) {
        codonCounts[codon] = 0;
    });

    _.each(sequences, function(seqObj) {
        var codon = seqObj.sequence;
        if (isValidCodon(codon) && codonCounts.hasOwnProperty(codon)) {
            codonCounts[codon]++;
        }
    });

    return codonCounts;
}

// Utility: Check if codon is valid
function isValidCodon(codon) {
    return codon.length === 3 && /^[ACGT]{3}$/.test(codon);
}

// Utility: Return expected codons and degeneracy for an amino acid
function getExpectedCodons(aminoAcid) {
    var codonData = {
        "F": { codons: ["TTT", "TTC"], sets: [["TTT", "TTC"]] },
        "L": { codons: ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], sets: [["TTA", "TTG"], ["CTT", "CTC", "CTA", "CTG"]] },
        "I": { codons: ["ATT", "ATC", "ATA"], sets: [["ATT", "ATC", "ATA"]] },
        "M": { codons: ["ATG"], sets: [["ATG"]] },
        "V": { codons: ["GTT", "GTC", "GTA", "GTG"], sets: [["GTT", "GTC", "GTA", "GTG"]] },
        "S": { codons: ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], sets: [["TCT", "TCC", "TCA", "TCG"], ["AGT", "AGC"]] },
        "P": { codons: ["CCT", "CCC", "CCA", "CCG"], sets: [["CCT", "CCC", "CCA", "CCG"]] },
        "T": { codons: ["ACT", "ACC", "ACA", "ACG"], sets: [["ACT", "ACC", "ACA", "ACG"]] },
        "A": { codons: ["GCT", "GCC", "GCA", "GCG"], sets: [["GCT", "GCC", "GCA", "GCG"]] },
        "Y": { codons: ["TAT", "TAC"], sets: [["TAT", "TAC"]] },
        "H": { codons: ["CAT", "CAC"], sets: [["CAT", "CAC"]] },
        "Q": { codons: ["CAA", "CAG"], sets: [["CAA", "CAG"]] },
        "N": { codons: ["AAT", "AAC"], sets: [["AAT", "AAC"]] },
        "K": { codons: ["AAA", "AAG"], sets: [["AAA", "AAG"]] },
        "D": { codons: ["GAT", "GAC"], sets: [["GAT", "GAC"]] },
        "E": { codons: ["GAA", "GAG"], sets: [["GAA", "GAG"]] },
        "C": { codons: ["TGT", "TGC"], sets: [["TGT", "TGC"]] },
        "W": { codons: ["TGG"], sets: [["TGG"]] },
        "R": { codons: ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], sets: [["CGT", "CGC", "CGA", "CGG"], ["AGA", "AGG"]] },
        "G": { codons: ["GGT", "GGC", "GGA", "GGG"], sets: [["GGT", "GGC", "GGA", "GGG"]] }
    };

    var entry = codonData[aminoAcid];
    if (entry) {
        return {
            codons: entry.codons,
            degeneracy: entry.codons.length,
            sets: entry.sets
        };
    } else {
        return null;
    }
}

// Utility: Identify atypical codons based on position and bias type
function isAtypical(codon, aminoAcid, biasType) {

    var preferredCodonsByBias = {
        "A-rich": {
            "A": "GCA", "V": "GTA", "L": "TTA", "P": "CCA", "T": "ACA", "S": "TCA",
            "R": "AGA", "G": "GGA", "N": "AAT", "K": "AAA", "D": "GAT", "E": "GAA",
            "Q": "CAA", "H": "CAT", "F": "TTT", "Y": "TAT", "C": "TGT", "W": "TGG",
            "M": "ATG", "I": "ATA"
        },
        "T-rich": {
            "A": "GCT", "V": "GTT", "L": "CTT", "P": "CCT", "T": "ACT", "S": "TCT",
            "R": "CGT", "G": "GGT", "N": "AAT", "K": "AAA", "D": "GAT", "E": "GAA",
            "Q": "CAA", "H": "CAT", "F": "TTT", "Y": "TAT", "C": "TGT", "W": "TGG",
            "M": "ATG", "I": "ATT"
        },
        "C-rich": {
            "A": "GCC", "V": "GTC", "L": "CTC", "P": "CCC", "T": "ACC", "S": "TCC",
            "R": "CGC", "G": "GGC", "N": "AAC", "K": "AAG", "D": "GAC", "E": "GAG",
            "Q": "CAG", "H": "CAC", "F": "TTC", "Y": "TAC", "C": "TGC", "W": "TGG",
            "M": "ATG", "I": "ATC"
        },
        "G-rich": {
            "A": "GCG", "V": "GTG", "L": "CTG", "P": "CCG", "T": "ACG", "S": "TCG",
            "R": "CGG", "G": "GGG", "N": "AAC", "K": "AAG", "D": "GAC", "E": "GAG",
            "Q": "CAG", "H": "CAC", "F": "TTC", "Y": "TAC", "C": "TGC", "W": "TGG",
            "M": "ATG", "I": "ATC"
        }
    };

    var biasTable = preferredCodonsByBias[biasType];
    if (biasTable) {
        var preferred = biasTable[aminoAcid];
        if (preferred) {
            return codon !== preferred;
        }
    }

    // If biasType is not recognized or no preferred codon is defined, treat as not atypical
    return false;
}
