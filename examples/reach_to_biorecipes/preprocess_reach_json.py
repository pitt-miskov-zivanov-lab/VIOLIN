# Created by Difei May 2023

# Usage
###############
# python preprocess_reach_json.py
###############

# Function
###############
# Take REACH outputs(multiple JSON files in "json" folder), assemble them into a single and large JSON file
###############

import json
import os
import glob

reach_input_dir = 'json'
path = os.path.join(reach_input_dir, '*.txt')
files = glob.glob(path)

data_events = dict()
data_entities = dict()
data_sentences = dict()
for fn in files:
    pathname, _ = os.path.splitext(fn)
    basename = os.path.basename(pathname)

    # Opening JSON file
    events_fn = os.path.join(reach_input_dir, basename + '.uaz.events.json')
    entities_fn = os.path.join(reach_input_dir, basename + '.uaz.entities.json')
    sentences_fn = os.path.join(reach_input_dir, basename + '.uaz.sentences.json')

    with open(events_fn, 'r', encoding='utf-8') as json_file:
        events = json.load(json_file)
        for i, value in events.items():
            if i == 'frames':
                if i not in data_events:
                    data_events[i] = value
                else:
                    data_events[i].extend(value)

    with open(entities_fn, 'r', encoding='utf-8') as json_file:
        entities = json.load(json_file)
        for i, value in entities.items():
            if i == 'frames':
                if i not in data_entities:
                    data_entities[i] = value
                else:
                    data_entities[i].extend(value)

    with open(sentences_fn, 'r', encoding='utf-8') as json_file:
        sentences = json.load(json_file)
        for i, value in sentences.items():
            if i == 'frames':
                if i not in data_sentences:
                    data_sentences[i] = value
                else:
                    data_sentences[i].extend(value)


    data = dict()

    data["events"] = data_events
    data["entities"] = data_entities
    data["sentences"] = data_sentences

    with open('reach_output_summary.json', 'w') as output_file:
        json.dump(data, output_file, indent=4)
