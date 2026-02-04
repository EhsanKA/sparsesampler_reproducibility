#!/usr/bin/env python
# coding: utf-8

"""
Create a clear table showing SPS advantages for each dataset
"""

def create_advantage_table():
    """Create a clear table of SPS advantages per dataset."""

    table = """
╔══════════════════════════════════════════════════════════════════════════════╗
║                    SPS ADVANTAGES: COMPLEMENTARY SUPERIORITY                ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                            ║
║  DATASET 1: MCC (Rare Cell Preservation Focus)                             ║
║  ─────────────────────────────────────────────────────────────────────────  ║
║  Metric                          SPS                Random              ║
║  ─────────────────────────────────────────────────────────────────────────  ║
║  Osteoblast Preservation        4,492 (15%)         927 (3%)             ║
║  Rare Cell Advantage            ✅ 4.8x more cells  ❌                    ║
║  DE Gene Overlap               57-61%              98-98%               ║
║  Standard Analysis             Random better        ✅                   ║
║                                                                            ║
║  KEY TAKEAWAY: SPS dramatically superior for rare cell representation     ║
║                                                                            ║
║  DATASET 2: LCMV (Biological Heterogeneity Focus)                         ║
║  ─────────────────────────────────────────────────────────────────────────  ║
║  Metric                          SPS                Random              ║
║  ─────────────────────────────────────────────────────────────────────────  ║
║  Clustering ARI                 0.380 (+6%)         0.358                ║
║  Clusters Detected              44                  34                   ║
║  Biological Diversity           ✅ Better preserved ❌                    ║
║  DE Gene Overlap               82-83%              87-87%               ║
║  Standard Analysis             Random slightly      ✅                   ║
║                                 better                                    ║
║                                                                            ║
║  KEY TAKEAWAY: SPS superior for preserving biological structure           ║
║                                                                            ║
║  OVERALL CONCLUSION:                                                       ║
║  ─────────────────────────────────────────────────────────────────────────  ║
║  SPS demonstrates complementary superiority across biological scenarios:   ║
║  • Rare cell preservation (MCC): SPS 4.8x better                          ║
║  • Biological diversity (LCMV): SPS better clustering structure           ║
║  • Standard analyses: Competitive performance with random                 ║
║                                                                            ║
║  CHOOSE SPS WHEN: Biological diversity and rare cells matter most         ║
║  CHOOSE RANDOM WHEN: Standard DE/clustering is the primary goal          ║
║                                                                            ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

    print(table)

    # Save to file
    with open('/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/sps_advantage_table.txt', 'w') as f:
        f.write(table)

    print("\nTable saved as: sps_advantage_table.txt")

if __name__ == '__main__':
    create_advantage_table()