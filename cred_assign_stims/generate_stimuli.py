"""
This is code to generate stimuli of the desired length

For ophys type sessions, unexpected sequences or violations are included. These are 
every 30-90 seconds, and last 2-4 seconds for the visual flow, and 3-6 seconds for the Gabors.
For hab type sessions, no unexpected sequences or violations occur.

Everything is randomized for each session for each animal (ordering of stim 
types (gabors vs visual flow), ordering of visual flow directions (left, right), 
positions, sizes and orientations of Gabors, location of visual flow squares, and time
and duration of unexpected sequences or violations if these occur.
"""

import copy
import logging
import os
import random
import sys
import pickle
import glob as glob
from win32api import GetSystemMetrics

import numpy as np
from psychopy import monitors

from camstim import Window, Warp

import stimulus_params
from stimulus_params import MOVIE_PARAMS
from cred_assign_stims import SweepStimModif, unique_directory

# Configuration settings used in the Credit Assignment project
WIDTH = 52.0
DISTANCE = 15.0
SIZEPIX = [1920, 1200]

def get_cred_assign_monitor():
    # Monitor sizing specs used in the Credit Assignment project
    monitor = monitors.Monitor("CredAssign")
    monitor.setWidth(WIDTH)
    monitor.setDistance(DISTANCE)
    monitor.setSizePix(SIZEPIX)

    return monitor


def check_reproduce(monitor, fullscreen=False, raise_error=False):
    
    orig_widpix, orig_heipix = SIZEPIX
    orig_wid = WIDTH
    orig_dist = DISTANCE

    curr_widpix, curr_heipix = monitor.getSizePix()
    fsc_str = ""
    if fullscreen:
        curr_widpix, curr_heipix = [GetSystemMetrics(i) for i in [0, 1]]
        fsc_str = " (fullscreen)"
    curr_wid = monitor.getWidth()
    curr_dist = monitor.getDistance()

    names = ["width (pix)", "height (pix)", "width (cm)", "distance (cm)"]
    values = [(orig_widpix, curr_widpix), (orig_heipix, curr_heipix), 
        (orig_wid, curr_wid), (orig_dist, curr_dist)]

    accumulate_names = []
    for name, (orig, curr) in zip(names, values):
        if orig != curr:
            accumulate_names.append(name)

    if len(accumulate_names) != 0:
        verb = "does" if len(accumulate_names) == 1 else "do"
        msg = ("Current {}{} {} not match original used in experiment. "
            "Seeds will not allow exact Credit Assignment stimulus parameter "
            "reproduction.").format(", ".join(accumulate_names), fsc_str, verb)
        if raise_error:
            raise ValueError(msg)
        else:
            logging.warning(msg)


def generate_stimuli(session_params, seed=None, save_frames="", save_directory=".", 
                     monitor=None, fullscreen=False, warp=False, save_from_frame=0):
    """
    generate_stimuli(session_params)

    Required args:
        - session_params (dict): see run_generate_stimuli.SESSION_PARAMS_OPHYS for 
                                 required keys and description.
    
    Optional args:
        - seed (int)           : seed to use to initialize Random Number Generator. 
                                 If None, will be set randomly.
                                 default: None
        - save_frames (str)    : extension to use for saving frames (frames not saved 
                                 if "")
                                 default: ""
        - save_directory (str) : main directory in which to save frames
                                 default: "."
        - monitor (Monitor)    : Psychopy Monitor. If None, a default test Monitor is 
                                 initialized instead.
                                 default: None
        - fullscreen (bool)    : If True, overrides monitor size
                                 default: False
        - warp (bool)          : If True, image is warped
                                 default: False
        - save_from_frame (int): Frame as of which to start saving frames, if saving
                                 default: 0
    """

    #recapitulate = raw_input("Are you recapitulating a specific session? y/n: ")
    recapitulate = 'y'
    # Get Pickle from experimental run
     # Get Pickle from experimental run
    path = "C:/Users/Henry Denny/camstim/output/"
    #filename = raw_input("Which file are you recapitulating? ")
    #if filename == 'test':
    filename = '221019212612-full_pipeline_script'
    fnm_glob = path + filename + "*.pkl"

    fnm = glob.glob(fnm_glob)[0]

    with open(fnm, 'rb') as file:
        stim_data = pickle.load(file)

    print(filename)

    # Setup session_structure dictionary, into which we will import the structure of the
    # experimental run
    session_structure = {
     #   "seed": "None",
        'size': 0,
        'display_sequence': {},
        'grt_sweep_order': [],
        'grt_sweep_table': [],
      #  "prev_session_params": 0,
    }

    # Get previous session params
    for i in range(len(stim_data['stimuli'])):
        if 'stim_params' in stim_data['stimuli'][i]:
            session_params = stim_data['stimuli'][i]['stim_params']['session_params']
            retrieved_rng = stim_data['stimuli'][i]['stim_params']['session_params']['rng']
            print(retrieved_rng.get_state()[1][1])
            session_structure['size'] = stim_data['monitor']['sizepix']

    # Get seed
    session_structure['seed'] = stim_data['stimuli'][6]['stim_params']['session_params']['seed']

    # Get movie display_sequences
    keys = []
    habkeys = ['m00', 'm10', 'm20']
    ophyskeys = ['m00', 'm01', 'm02','m03', 'm10', 'm11', 'm12', 'm13', 'm20', 'm21', 'm22', 'm23']
    movcount = 0
    varicount = 0

    # Find which stimuli contain movies
    for i in np.arange(len(stim_data['stimuli'])):
        if 'movie_path' in stim_data['stimuli'][i]:
            keys.append(i)
    
    # Retrieve movie clip display times.    
    if session_params['type'] == 'ophys':
        session_structure['display_sequence'] = dict.fromkeys(ophyskeys)
        for i, j in enumerate(keys):
            session_structure['display_sequence'][str(ophyskeys[i])] = stim_data['stimuli'][j]['display_sequence']
    elif session_params['type'] == 'hab':
        session_structure['display_sequence'] = dict.fromkeys(habkeys)
        for i, j in enumerate(keys):
            session_structure['display_sequence'][str(habkeys[i])] = stim_data['stimuli'][i]['display_sequence']
        
    #Get gratings parameters
    session_structure['grt_sweep_order'] = stim_data['stimuli'][len(stim_data['stimuli'])-1]['sweep_order']
    session_structure['grt_sweep_table'] = stim_data['stimuli'][len(stim_data['stimuli'])-1]['sweep_table']

    # Record orientations of gabors at each sweep (LEAVE AS TRUE)
    recordOris = True

    # Record positions of squares at all times (LEAVE AS TRUE)
    recordPos = True
            
    # create a monitor
    if monitor is None:
        get_cred_assign_monitor()
    
    check_reproduce(monitor, fullscreen=fullscreen, raise_error=False)

    spare = session_params["rng"].get_state()
    
    #if recapitulate == 'n':
            # randomly set a seed for the session
        #session_params["seed"] = random.randint(1, 10000)
    #else:
    #    session_params["seed"] = session_structure['seed']
    logging.info("Seed: {}".format(session_params["seed"]))
    #session_params["rng"] = np.random.RandomState(session_params["seed"])
    # recapRS = np.random.RandomState(session_params["seed"])
    # session_params["rng"] = retrieved_rng
    # print(retrieved_rng.get_state()[1][1])

    session_params['rng'] = np.random.RandomState(session_params['seed'])
    # print('sum rs = '+str(np.subtract(session_params["rng"].get_state()[1], recapRS.get_state()[1][1])))

    # print('Recap RandomState = ' + str(session_params["rng"].get_state()[1][1]))
    # print(spare[1][1])
    # print(np.subtract(session_params["rng"].get_state()[1], spare[1]))
    # print(np.subtract(recapRS.get_state()[1][1], spare[1]))
    # print(recapRS.get_state()[1][1])



    # check session params add up to correct total time
    tot_calc = session_params["pre_blank"] + session_params["post_blank"] + \
               2 * session_params["inter_blank"] + session_params["gab_dur"] + \
               2 * session_params["sq_dur"]
    if tot_calc != session_params["session_dur"]:
        logging.warning("Session expected to add up to {} s, but adds up to {} s."
              .format(session_params["session_dur"], tot_calc))

    # Create display window
    if recapitulate == 'y':
        window_kwargs = {
            "fullscr": fullscreen,
            "size"   : [stim_data['wwidth'], stim_data['wheight']], # May return an error due to size. Ignore.
            "monitor": monitor, # Will be set to a gamma calibrated profile by MPE
            "screen" : 0,
            }
    else: 
        window_kwargs = {
            "fullscr": fullscreen,
            "size"   : monitor.getSizePix(), # May return an error due to size. Ignore.
            "monitor": monitor, # Will be set to a gamma calibrated profile by MPE
            "screen" : 0,
        }
    
    if warp:
        window_kwargs["warp"] = Warp.Spherical

    if recapitulate == 'y':
        window_kwargs["warp"] = Warp.Spherical
        window_kwargs["fullscr"] = False     
        #window_kwargs["size"] = session_structure['size']

    window = Window(**window_kwargs)

    stim_order = []
    sq_order = []
    gab_order = []
    rot_gab_order = []
    mov_order = []
    grt_order = []
    gab_block_order = []
    interpreter = 0
   
    # initialize the stimuli
    if session_params['gab_dur'] != 0:
        gb_1 = stimulus_params.init_run_gabors(window, session_params.copy(), recordOris, surp=2)
        
        # share positions and sizes
        gb_2_session_params = session_params.copy()
        gb_2_session_params['possize'] = gb_1.stim_params['session_params']['possize']
        gb_2 = stimulus_params.init_run_gabors(window, gb_2_session_params, recordOris, surp=2)
        
        stim_order.append('g')
        gab_order = [1, 2]
        gab_block_order = [1, 2]
    if session_params['rot_gab_dur'] != 0:
        rgb_1 = stimulus_params.init_rotate_gabors(window, gb_2_session_params, recordOris, surp=2)
        
        # share positions and sizes from original Gabors. Keeps possize the same
        rgb_2 = stimulus_params.init_rotate_gabors(window, gb_2_session_params, recordOris, surp=2)
        rot_gab_order = [1, 2]
    if session_params['sq_dur'] != 0:
        sq_left = stimulus_params.init_run_squares(window, 'left', session_params.copy(), recordPos)
        sq_right = stimulus_params.init_run_squares(window, 'right', session_params.copy(), recordPos)
        stim_order.append('b')
        sq_order = ['l', 'r']
    if session_params['movie_dur'] != 0:
        mov, propblocks = stimulus_params.init_run_movies(window, session_params.copy(), MOVIE_PARAMS, 1, session_params['movie_folder'])
        stim_order.append('m')
        mov_order = np.arange(MOVIE_PARAMS['movie_n'])
    if session_params['gratings_dur'] != 0:
        grt = stimulus_params.init_run_gratings(window, session_params.copy())
        stim_order.append('grt')

    print(stim_order)

    # initialize display order and times
    session_params['rng'].shuffle(stim_order) # in place shuffling
    session_params['rng'].shuffle(sq_order) # in place shuffling
    session_params['rng'].shuffle(gab_order) # in place shuffling
    session_params['rng'].shuffle(rot_gab_order) # in place shuffling
    session_params['rng'].shuffle(gab_block_order) # in place shuffling

    print(stim_order)
    
    displayorder = {}
    if session_params['type'] == 'ophys':    
        for i in np.arange(MOVIE_PARAMS['vids_per_block']):
            displayorder[str(i)] = []
    elif session_params['type'] == 'hab':
        for i in np.arange(0, MOVIE_PARAMS['vids_per_block'], 4):
            displayorder[str(i)] = []

    start = session_params["pre_blank"] # initial blank
    stimuli = []
    ID = 0
    
    for i in stim_order:
        if i == 'g':
            for l in gab_block_order:
                if l == 1:
                    for j in gab_order:
                        if j == 1:
                            stimuli.append(gb_1)
                            gb_1.set_display_sequence([(start, start+session_params['gab_dur'])])
                        elif j == 2:
                            stimuli.append(gb_2)
                            gb_2.set_display_sequence([(start, start+session_params['gab_dur'])])
                        # update the new starting point for the next stim
                        start += session_params['gab_dur'] 
                        print('t after fg' + str(start/60))
                elif l == 2:
                    for j in rot_gab_order:
                        if j == 1:
                            stimuli.append(rgb_1)
                            rgb_1.set_display_sequence([(start, start+session_params['rot_gab_dur'])])
                        elif j == 2:
                            stimuli.append(rgb_2)
                            rgb_2.set_display_sequence([(start, start+session_params['rot_gab_dur'])])
                        start += session_params['gab_dur']
                        print('t after rg' + str(start/60))
                        # update the new starting point for the next stim
                start += session_params['inter_blank']
        elif i == 'b':
            for j in sq_order:
                if j == 'l':
                    stimuli.append(sq_left)
                    sq_left.set_display_sequence([(start, start+session_params['sq_dur'])])
                elif j == 'r':
                    stimuli.append(sq_right)
                    sq_right.set_display_sequence([(start, start+session_params['sq_dur'])])
                # update the new starting point for the next stim
                start += session_params['sq_dur'] + session_params['inter_blank'] 
        elif i == 'm':
            if recapitulate == 'n':    
                if session_params['type'] == 'ophys':
                    for ii in np.arange(session_params['movie_blocks']):
                        propblocksshuf = np.random.permutation(propblocks)
                        for j in propblocksshuf:
                            displayorder[str(j)].append((start, start+(MOVIE_PARAMS['movie_len'])-1))
                            start += MOVIE_PARAMS['movie_len']
                    for j in np.arange(MOVIE_PARAMS['vids_per_block']):
                        mov[str(j)].set_display_sequence(displayorder[str(j)])
                        stimuli.append(mov[str(j)])
                elif session_params['type'] == 'hab':
                    for ii in np.arange(session_params['movie_blocks']):
                        propblocksshuf = np.random.permutation(propblocks)
                        for j in propblocksshuf:
                            displayorder[str(j)].append((start, start+(MOVIE_PARAMS['movie_len'])-1))
                            start += MOVIE_PARAMS['movie_len']
                    for j in np.arange(0, MOVIE_PARAMS['vids_per_block'], 4):
                        mov[str(j)].set_display_sequence(displayorder[str(j)])
                        stimuli.append(mov[str(j)])
            elif recapitulate == 'y':
                if session_params['type'] == 'ophys':
                    for l, j in enumerate(keys):
                        mov[str(l)].set_display_sequence(session_structure['display_sequence'][ophyskeys[l]])
                        stimuli.append(mov[str(l)])
                if session_params['type'] == 'hab':
                    for l, j in enumerate(habkeys):
                        print(l)
                        print(habkeys[l])
                        print(ID)
                        mov[str(ID)].set_display_sequence(session_structure['display_sequence'][habkeys[l]])
                        stimuli.append(mov[str(ID)])
                        ID = ID + 4 
            start += MOVIE_PARAMS['movie_len']*MOVIE_PARAMS['vids_per_block']*session_params['movie_blocks']
            start += session_params['inter_blank']
            # update the new starting point for the next stim
    if session_params['gratings_dur'] != 0:
        if recapitulate == 'y':
            grt.sweep_order = session_structure['grt_sweep_order']
            grt.sweep_table = session_structure['grt_sweep_table']
        grt.set_display_sequence([(start, (start + session_params['gratings_dur']*14.7))])
        stimuli.append(grt)

    # prepare path for file saving
    frames_path = ""
    if save_frames:
        frames_directory = os.path.join(
            save_directory, "frames_{}".format(str(session_params["seed"])))
        frames_path = os.path.join(frames_directory, "frame_.{}".format(save_frames))

    ss = SweepStimModif(
        window=window,
        stimuli=stimuli,
        post_blank_sec=session_params["post_blank"],
        params={},  # will be set by MPE to work on the rig
        frames_output=frames_path,
        save_from_frame=save_from_frame,
        name=session_params["seed"],
        warp=warp,
        set_brightness=False # skip setting brightness
        )

    # catch system exit 0 thrown by ss._finalize()
    try:
        ss.run()
    except SystemExit as cm:
        if cm.code != 0:
            sys.exit(cm.code)
        logging.warning("Ignoring automatic post-sweep system exit.")


