#!/usr/bin/env python3
"""
Generate a synthetic test ntuple for validating PPSG4NtupleToRecHitProducer.

This script creates a ROOT file with the same structure as the G4 simulation output,
with known hit positions that can be used to verify the reconstruction chain.

Test cases:
1. Straight-line protons through detector center -> tracks with zero slope
2. Protons with known transverse slopes -> fitted tracks should recover the slope
3. Events with missing planes -> tests hit efficiency handling
"""

import ROOT
import math
import argparse

# Detector parameters (matching G4 simulation)
NUM_PLANES = 6
PIXEL_PITCH_X = 0.1  # mm (100 um)
PIXEL_PITCH_Y = 0.15  # mm (150 um)
NUM_ROWS = 160
NUM_COLS = 104
ACTIVE_HALF_X = (NUM_ROWS * PIXEL_PITCH_X) / 2.0  # 8.0 mm
ACTIVE_HALF_Y = (NUM_COLS * PIXEL_PITCH_Y) / 2.0  # 7.8 mm

# Plane z-positions (relative, in mm) - typical spacing
PLANE_Z = [0.0, 4.0, 8.0, 12.0, 16.0, 20.0]


def position_to_pixel(x, y):
    """Convert local position (mm) to pixel indices."""
    row = int((x + ACTIVE_HALF_X) / PIXEL_PITCH_X)
    col = int((y + ACTIVE_HALF_Y) / PIXEL_PITCH_Y)
    row = max(0, min(row, NUM_ROWS - 1))
    col = max(0, min(col, NUM_COLS - 1))
    return row, col


def create_test_ntuple(output_file, num_events=100):
    """Create a test ntuple with various track configurations."""
    f = ROOT.TFile(output_file, 'RECREATE')
    f.SetCompressionLevel(0)

    # Create scoring plane tree (truth info)
    truth_tree = ROOT.TTree('plane0', 'Scoring plane 0 hits')
    truth_x = ROOT.std.vector('double')()
    truth_y = ROOT.std.vector('double')()

    # Use arrays for simple branch binding
    t_x = [0.0]
    t_y = [0.0]
    t_theta_x = [0.0]
    t_theta_y = [0.0]
    t_p = [0.0]
    t_energy = [0.0]
    t_eventID = [0]

    truth_tree.Branch('x', t_x, 'x/D')
    truth_tree.Branch('y', t_y, 'y/D')
    truth_tree.Branch('theta_x', t_theta_x, 'theta_x/D')
    truth_tree.Branch('theta_y', t_theta_y, 'theta_y/D')
    truth_tree.Branch('p', t_p, 'p/D')
    truth_tree.Branch('primary_energy', t_energy, 'primary_energy/D')
    truth_tree.Branch('eventID', t_eventID, 'eventID/I')

    # Create rechit trees (one per plane)
    rechit_trees = []
    rechit_data = []
    for i in range(NUM_PLANES):
        tree = ROOT.TTree(f'rechits_plane{i+1}', f'RecHits for pixel plane {i+1}')
        data = {
            'eventID': [0],
            'row': [0],
            'col': [0],
            'x': [0.0],
            'y': [0.0],
            'edep': [0.0],
            'nhits': [0]
        }
        tree.Branch('eventID', data['eventID'], 'eventID/I')
        tree.Branch('row', data['row'], 'row/I')
        tree.Branch('col', data['col'], 'col/I')
        tree.Branch('x', data['x'], 'x/D')
        tree.Branch('y', data['y'], 'y/D')
        tree.Branch('edep', data['edep'], 'edep/D')
        tree.Branch('nhits', data['nhits'], 'nhits/I')
        rechit_trees.append(tree)
        rechit_data.append(data)

    # Generate test events
    for evt in range(num_events):
        # Vary track parameters based on event number
        if evt < num_events // 4:
            # Type 1: Straight through center (zero slope)
            x0, y0 = 0.0, 0.0
            slope_x, slope_y = 0.0, 0.0
        elif evt < num_events // 2:
            # Type 2: Small positive slope
            x0, y0 = -1.0, -0.5
            slope_x = 0.01  # 10 mrad
            slope_y = 0.005  # 5 mrad
        elif evt < 3 * num_events // 4:
            # Type 3: Small negative slope
            x0, y0 = 1.0, 0.5
            slope_x = -0.008
            slope_y = -0.003
        else:
            # Type 4: Varying parameters
            x0 = -2.0 + 4.0 * (evt - 3 * num_events // 4) / (num_events // 4)
            y0 = -1.0 + 2.0 * (evt - 3 * num_events // 4) / (num_events // 4)
            slope_x = 0.002 * math.sin(evt * 0.1)
            slope_y = 0.001 * math.cos(evt * 0.15)

        # Momentum (typical forward proton: ~6500 GeV)
        momentum = 6500.0 - 100.0 * (evt % 10)  # Slight variation

        # Fill truth tree
        t_x[0] = x0
        t_y[0] = y0
        t_theta_x[0] = slope_x * 1000.0  # Convert to mrad
        t_theta_y[0] = slope_y * 1000.0
        t_p[0] = momentum
        t_energy[0] = math.sqrt(momentum**2 + 0.938**2)
        t_eventID[0] = evt
        truth_tree.Fill()

        # Skip some planes for some events (to test missing hits)
        skip_planes = []
        if evt % 20 == 19:
            skip_planes = [2]  # Skip plane 3
        if evt % 30 == 29:
            skip_planes = [0, 5]  # Skip planes 1 and 6

        # Fill rechit trees for each plane
        for plane_idx in range(NUM_PLANES):
            if plane_idx in skip_planes:
                continue

            # Calculate hit position at this plane
            z = PLANE_Z[plane_idx]
            hit_x = x0 + slope_x * z
            hit_y = y0 + slope_y * z

            # Check if within active area
            if abs(hit_x) > ACTIVE_HALF_X or abs(hit_y) > ACTIVE_HALF_Y:
                continue

            # Convert to pixel indices
            row, col = position_to_pixel(hit_x, hit_y)

            # Fill rechit tree
            data = rechit_data[plane_idx]
            data['eventID'][0] = evt
            data['row'][0] = row
            data['col'][0] = col
            data['x'][0] = hit_x
            data['y'][0] = hit_y
            data['edep'][0] = 100.0  # Arbitrary energy deposit in keV
            data['nhits'][0] = 1
            rechit_trees[plane_idx].Fill()

    # Write and close
    f.cd()
    truth_tree.Write()
    for tree in rechit_trees:
        tree.Write()
    f.Close()

    print(f'Created test ntuple: {output_file}')
    print(f'  Events: {num_events}')
    print(f'  Scoring plane entries: {truth_tree.GetEntries()}')
    for i, tree in enumerate(rechit_trees):
        print(f'  rechits_plane{i+1} entries: {tree.GetEntries()}')


def main():
    parser = argparse.ArgumentParser(description='Generate test ntuple for PPS G4 bridge')
    parser.add_argument('-o', '--output', default='test_ntuple.root',
                        help='Output ROOT file path')
    parser.add_argument('-n', '--nevents', type=int, default=100,
                        help='Number of events to generate')
    args = parser.parse_args()

    create_test_ntuple(args.output, args.nevents)


if __name__ == '__main__':
    main()
