import ROOT
import argparse
import os
import logging
import tkinter as tk
from itertools import product

ARMS = [0,1]
STATIONS = [0,2] 
PLANES = [0,1,2,3,4,5] 

def makeCanvasWithAspectRatio(name, title, w, h):

    root = tk.Tk()
    root.withdraw()  # Hide the main window
    displayW = root.winfo_screenwidth()
    displayH = root.winfo_screenheight()
    finalW = displayW - 150
    finalH = finalW * h / w

    h_limited = finalH > (displayH - 15)
    if h_limited:
        finalH = displayH - 150
        finalW = finalH * w / h

    return ROOT.TCanvas(name, title, 100, 100, int(finalW), int(finalH))

def plotPlaneEfficiencies(file, canvases, hists):
    previous_station = -1
    for arm, station, plane in product(ARMS, STATIONS, PLANES):
        logging.debug(f'Processing Arm {arm}, Station {station}, Plane {plane}')

        if station != previous_station:
                        # Create a new canvas for the new station
            canvas = makeCanvasWithAspectRatio(f'c_plane_arm{arm}_station{station}', f'Plane efficiencies - Arm {arm} Station {station}', 14, 9)
            canvases.append(canvas)
            canvas.Divide(3, 2)

        canvas.cd(plane+1)
        canvas.GetPad(plane+1).SetRightMargin(0.15)
        hist_name = f'DQMData/Run 999999/Arm{arm}/Run summary/st{station}/rp3/' \
                    f'h2EfficiencyMap_arm{arm}_st{station}_rp3_pl{plane}'
        hist = file.Get(hist_name).Clone()
        if not hist:
            logging.warning(f'Histogram {hist_name} not found in file {file.GetName()}. Skipping.')
            continue
        hists.append(hist)
        sector = '45' if arm == 0 else '56'
        station_name = '210' if station == 0 else '220'
        hist.SetTitle(f'{sector}-{station_name} plane {plane}')
        hist.GetXaxis().SetTitle('x (mm)')
        hist.GetXaxis().SetTitleSize(0.04)
        hist.GetYaxis().SetTitle('y (mm)')
        hist.GetYaxis().SetTitleSize(0.04)
        hist.GetZaxis().SetRangeUser(0, 1)
        hist.GetZaxis().SetTitleSize(0.04)
        hist.GetZaxis().SetTitle('Efficiency')
        hist.GetZaxis().SetTitleOffset(1.2)
        hist.Draw('COLZ')
        canvas.SetLogz(True)
        canvas.Update()

        previous_station = station

    return canvases, hists

def plot(efficiencyFile):
    # Check if the efficiency file exists
    if not os.path.isfile(efficiencyFile):
        logging.fatal(f'File {efficiencyFile} does not exist.')
        return

    ROOT.gStyle.SetOptStat(0)
    # Open the efficiency file
    file = ROOT.TFile(efficiencyFile, 'READ')
    if not file or file.IsZombie():
        logging.fatal(f'Could not open file {efficiencyFile}.')
        return

    # Lists to hold canvases and histograms
    canvases = []
    hists = []

    # File types:
    # - analysis: plane efficiency and multirp efficiency
    # - reference: unit efficiency -> not yet supported

    # Determine the type of the file
    # Check if the file contains the plane efficiency histograms
    file_type = None
    template_eff_hist_name = f'DQMData/Run 999999/Arm{ARMS[0]}/Run summary/st{STATIONS[0]}/rp3/' \
                            f'h2EfficiencyMap_arm{ARMS[0]}_st{STATIONS[0]}_rp3_pl{PLANES[0]}'
    if file.Get(template_eff_hist_name):
        logging.info(f'File {efficiencyFile} is an analysis file.')
        file_type = 'analysis'

    if file_type is None:
        logging.fatal(f'File {efficiencyFile} is not a valid efficiency file.')
        return
    elif file_type == 'analysis':
        logging.info(f'Processing analysis file {efficiencyFile}.')
        logging.warning('Only plane efficiency histograms are supported for plotting at the moment.')
        canvases, hists = plotPlaneEfficiencies(file, canvases, hists)
    
    return canvases, hists

def parseCommand():
    # Create the top-level command parser
    parser=argparse.ArgumentParser(
        description='This tool helps with handling the PPS pixel mask files.',
        epilog='Made by A. Bellora, 10/11/2022'
    )
    parser.add_argument('file', type=str, help='Path to the efficiency file')
    parser.add_argument('--debug', action='store_true', default=False, help='Enable debug mode (default: False)')

    return parser.parse_args()

if __name__ == '__main__':
    args = parseCommand()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG,
                            format='[%(asctime)s - '+os.path.basename(__file__)+' - %(levelname)s] %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
        
        logging.info('Debug mode enabled')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='[%(asctime)s - '+os.path.basename(__file__)+' - %(levelname)s] %(message)s',datefmt='%Y-%m-%d %H:%M:%S')

    canvases, hists = plot(args.file)
    input('Press Enter to exit...')