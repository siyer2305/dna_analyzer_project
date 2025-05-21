from flask import Flask, request, jsonify, render_template
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import re

app = Flask(__name__, template_folder='../templates', static_folder='../static')

def analyze_sequence(dna_sequence):
    # Clean the input sequence
    cleaned_sequence = re.sub(r'\s+', '', dna_sequence).upper()  # Remove whitespace and convert to uppercase
    
    # Input validation (more permissive)
    if not re.fullmatch(r'^[ATCGatcg]+$', cleaned_sequence):
        return jsonify({'error': 'Invalid DNA sequence. Only A, T, C, and G are allowed (case insensitive).'}), 400

    # Calculate GC content
    seq = Seq(cleaned_sequence)
    gc_content = gc_fraction(seq) * 100

    # Identify restriction sites
    restriction_enzymes = {
        "GAATTC": {"name": "EcoRI", "sequence": "GAATTC"},
        "CCGG": {"name": "HpaII", "sequence": "CCGG"},
        "GACGTC": {"name": "SalI", "sequence": "GACGTC"},
        "AAGCTT": {"name": "HindIII", "sequence": "AAGCTT"},
        "GGATCC": {"name": "BamHI", "sequence": "GGATCC"}
    }
    
    found_sites = {}
    for sequence, enzyme_info in restriction_enzymes.items():
        indices = [i + 1 for i in range(len(cleaned_sequence) - len(sequence) + 1) 
                  if cleaned_sequence[i:i + len(sequence)] == sequence]
        if indices:
            found_sites[enzyme_info["name"]] = {
                "sequence": enzyme_info["sequence"],
                "positions": indices
            }

    # ORF prediction function
    def find_orfs(sequence, strand, offset=0):
        orfs = []
        start_codon = "ATG"
        stop_codons = ["TAA", "TAG", "TGA"]
        min_orf_length = 60
        
        for frame in range(3):
            for i in range(frame, len(sequence) - 2, 3):
                if sequence[i:i+3] == start_codon:
                    for j in range(i + 3, len(sequence) - 2, 3):
                        codon = sequence[j:j+3]
                        if codon in stop_codons:
                            orf_sequence = sequence[i:j+3]
                            orf_length = len(orf_sequence)
                            if orf_length >= min_orf_length:
                                start = i + 1 + offset if strand == '+' else len(cleaned_sequence) - j - 2
                                end = j + 3 + offset if strand == '+' else len(cleaned_sequence) - i
                                orfs.append({
                                    "start": start,
                                    "end": end,
                                    "length": orf_length,
                                    "strand": strand,
                                    "frame": frame + 1,
                                    "sequence": orf_sequence
                                })
                            break
        return orfs

    # Find ORFs on both strands
    orfs = find_orfs(cleaned_sequence, '+')
    reverse_complement = str(seq.reverse_complement())
    orfs += find_orfs(reverse_complement, '-')

    return {
        'gc_content': gc_content,
        'restriction_sites': found_sites,
        'orfs': sorted(orfs, key=lambda x: x['length'], reverse=True),
        'dna_sequence': cleaned_sequence,
        'success': True
    }

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze():
    try:
        data = request.get_json()
        if not data or 'dna_sequence' not in data:
            return jsonify({'error': 'No DNA sequence provided'}), 400
            
        dna_sequence = data.get('dna_sequence', '').strip()
        if not dna_sequence:
            return jsonify({'error': 'Empty DNA sequence'}), 400
            
        result = analyze_sequence(dna_sequence)
        return jsonify(result)
        
    except Exception as e:
        return jsonify({'error': f'An error occurred: {str(e)}'}), 500

if __name__ == '__main__':
    app.run(debug=True)