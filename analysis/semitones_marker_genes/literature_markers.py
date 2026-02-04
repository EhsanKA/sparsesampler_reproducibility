#!/usr/bin/env python
# coding: utf-8

"""
Literature-based marker genes for MCC dataset cell types.

These marker genes are compiled from CellMarker 2.0 and PanglaoDB databases,
as well as published single-cell RNA-seq studies.

References:
- Hu C., et al. (2023). CellMarker 2.0. Nucleic Acids Research, 51(D1), D870–D876.
- Franzén O., et al. (2019). PanglaoDB. Database, 2019, baz046.
- Loh K.M., et al. (2016). Cell, 166(2), 451-467. (PMC5474394)
"""

# Literature-based marker genes for each cell type in the MCC dataset
LITERATURE_MARKERS = {
    'osteoblast': {
        'markers': [
            'RUNX2',   # Master transcription factor for osteoblast differentiation
            'SP7',     # Osterix - essential for bone formation
            'ALPL',    # Alkaline Phosphatase - early osteoblast marker
            'COL1A1',  # Collagen type I - major matrix protein
            'BGLAP',   # Osteocalcin - late marker, mature osteoblasts
            'SPARC',   # Osteonectin
            'SPP1',    # Osteopontin
            'IBSP',    # Bone sialoprotein
        ],
        'references': [
            'CellMarker 2.0 (Hu et al., 2023)',
            'PanglaoDB (Franzén et al., 2019)',
        ]
    },
    'chondrocyte': {
        'markers': [
            'SOX9',    # Key transcription factor for chondrogenic lineage
            'SOX5',    # Chondrogenic transcription factor
            'SOX6',    # Chondrogenic transcription factor
            'COL2A1',  # Collagen type II - cartilage matrix
            'ACAN',    # Aggrecan - major proteoglycan
            'COL10A1', # Hypertrophic chondrocyte marker
            'COL11A1', # Cartilage collagen
            'COMP',    # Cartilage oligomeric matrix protein
        ],
        'references': [
            'CellMarker 2.0 (Hu et al., 2023)',
            'Ji Q., et al. (2019). Ann Rheum Dis. (PMC7724342)',
        ]
    },
    'fibroblast': {
        'markers': [
            'VIM',     # Vimentin - mesenchymal cytoskeleton
            'PDGFRA',  # Growth factor receptor
            'PDGFRB',  # Growth factor receptor
            'COL1A1',  # Collagen type I
            'COL1A2',  # Collagen type I
            'THY1',    # CD90 - surface marker
            'FAP',     # Fibroblast activation protein
            'DCN',     # Decorin
            'LUM',     # Lumican
        ],
        'references': [
            'CellMarker 2.0 (Hu et al., 2023)',
            'PanglaoDB (Franzén et al., 2019)',
        ]
    },
    'mesodermal cell': {
        'markers': [
            'TBXT',    # Brachyury/T - classic mesoderm marker
            'MIXL1',   # Early mesendoderm marker
            'MESP1',   # Cardiac mesoderm specification
            'EOMES',   # Mesendoderm transcription factor
            'KDR',     # VEGFR2 - mesodermal precursors
            'GSC',     # Goosecoid
            'TBX6',    # T-box transcription factor
        ],
        'references': [
            'CellMarker 2.0 (Hu et al., 2023)',
            'PMID:19134196',
        ]
    },
    'lateral mesodermal cell': {
        'markers': [
            'HAND1',   # TF, LPM specific (~98-100% of LPM cells)
            'FOXF1',   # Central LPM marker (~98-100% of LPM cells)
            'GATA4',   # Downstream of BMP4 and FOXF1
            'BMP4',    # Key signaling molecule in LPM
            'WT1',     # Splanchnic mesenchyme marker
            'TBX18',   # Splanchnic mesenchyme marker
            'HAND2',   # Heart/LPM marker
        ],
        'references': [
            'Loh K.M., et al. (2016). Cell, 166(2), 451-467. (PMC5474394)',
            'Reactome Pathway: R-HSA-9758920',
        ]
    }
}

# Alternative names/aliases for marker genes (some datasets use different names)
GENE_ALIASES = {
    'TBXT': ['T', 'BRACHYURY'],
    'SP7': ['OSX', 'OSTERIX'],
    'BGLAP': ['OSTEOCALCIN', 'OCN'],
    'SPP1': ['OPN', 'OSTEOPONTIN'],
    'THY1': ['CD90'],
    'KDR': ['VEGFR2', 'FLK1'],
}


def get_all_markers():
    """Get a flat list of all unique marker genes."""
    all_markers = set()
    for cell_type, data in LITERATURE_MARKERS.items():
        all_markers.update(data['markers'])
    return sorted(list(all_markers))


def get_markers_for_celltype(cell_type):
    """Get marker genes for a specific cell type."""
    # Normalize cell type name
    cell_type_lower = cell_type.lower().strip()
    
    for ct, data in LITERATURE_MARKERS.items():
        if ct.lower() == cell_type_lower:
            return data['markers']
    
    # Try partial matching
    for ct, data in LITERATURE_MARKERS.items():
        if cell_type_lower in ct.lower() or ct.lower() in cell_type_lower:
            return data['markers']
    
    return []


def get_marker_celltype_mapping():
    """Get a dictionary mapping each marker to its cell type(s)."""
    marker_to_celltype = {}
    for cell_type, data in LITERATURE_MARKERS.items():
        for marker in data['markers']:
            if marker not in marker_to_celltype:
                marker_to_celltype[marker] = []
            marker_to_celltype[marker].append(cell_type)
    return marker_to_celltype


if __name__ == '__main__':
    print("Literature-based marker genes for MCC dataset:")
    print("=" * 60)
    for cell_type, data in LITERATURE_MARKERS.items():
        print(f"\n{cell_type.upper()}")
        print(f"  Markers: {', '.join(data['markers'])}")
        print(f"  References: {'; '.join(data['references'])}")
    
    print("\n" + "=" * 60)
    print(f"Total unique markers: {len(get_all_markers())}")
