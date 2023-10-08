import os
import matplotlib

import matplotlib.pyplot as plt
files = []
[files.append("data/PSD/mode_0/fm_demod/" + x) for x in os.listdir("data/PSD/mode_0/fm_demod")]
[files.append("data/PSD/mode_0/monodata/" + x) for x in os.listdir("data/PSD/mode_0/monodata")]
[files.append("data/PSD/mode_0/i_data/" + x) for x in os.listdir("data/PSD/mode_0/i_data")]
[files.append("data/PSD/mode_0/q_data/" + x) for x in os.listdir("data/PSD/mode_0/q_data")]


#[files.append("data/PSD/mode_1/fm_demod/" + x) for x in os.listdir("data/PSD/mode_1/fm_demod")]
#[files.append("data/PSD/mode_1/monodata/" + x) for x in os.listdir("data/PSD/mode_1/monodata")]
#[files.append("data/PSD/mode_1/i_data/" + x) for x in os.listdir("data/PSD/mode_1/i_data")]
#[files.append("data/PSD/mode_1/q_data/" + x) for x in os.listdir("data/PSD/mode_1/q_data")]

[files.append("data/PSD/mode_2/fm_demod/" + x) for x in os.listdir("data/PSD/mode_2/fm_demod")]
[files.append("data/PSD/mode_2/monodata/" + x) for x in os.listdir("data/PSD/mode_2/monodata")]
[files.append("data/PSD/mode_2/i_data/" + x) for x in os.listdir("data/PSD/mode_2/i_data")]
[files.append("data/PSD/mode_2/q_data/" + x) for x in os.listdir("data/PSD/mode_2/q_data")]

#[files.append("data/PSD/mode_3/fm_demod/" + x) for x in os.listdir("data/PSD/mode_3/fm_demod")]
#[files.append("data/PSD/mode_3/monodata/" + x) for x in os.listdir("data/PSD/mode_3/monodata")]
#[files.append("data/PSD/mode_3/i_data/" + x) for x in os.listdir("data/PSD/mode_3/i_data")]
#[files.append("data/PSD/mode_3/q_data/" + x) for x in os.listdir("data/PSD/mode_3/q_data")]


for file_name in files:
    #file_name = "src/monodata.txt"
    if file_name.endswith(".txt"):
        with open(file_name) as file:
            data = file.readlines()
            dataX = [float(x.split(", ")[0]) for x in data]
            dataY = [float(x.split(", ")[1]) for x in data]

        plt.plot(dataX, dataY)
        plt.savefig(file_name + ".png")
        plt.clf()
        print(f"plotted {file_name}")
    