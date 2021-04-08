import os
import pickle
import json


def write_dict_json(name, dir, d):
    if not os.path.exists(dir):
        os.mkdir(dir)

    path = os.path.join(dir, name + ".json")
    with open(path, 'w') as file:
        json.dump(d, file)


def write_data_synteny(smg, sml, indexes, i, name):
    if not os.path.exists("temp_" + name):
        os.mkdir("temp_" + name)
    with open("temp_" + name + "/thread_" + str(i + 1) +
              "_smg.temp", "wb") as file:
        pickle.dump(smg, file)
    with open("temp_" + name + "/thread_" + str(i + 1) +
              "_sml.temp", "wb") as file:
        pickle.dump(sml, file)
    with open("temp_" + name + "/thread_" + str(i + 1) +
              "_indexes.temp", "wb") as file:
        pickle.dump(indexes, file)
