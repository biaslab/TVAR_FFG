import matplotlib.pyplot as plt
from scipy.io import wavfile
plt.style.use("ggplot")


fs, y = wavfile.read("data/speech/clean.wav")
fs, x = wavfile.read("data/speech/noised.wav")
fs, z = wavfile.read("data/speech/filtered.wav")


fig, axs = plt.subplots(3)
axs[0].specgram(y, Fs=fs)
axs[0].set_title("Clean")
axs[1].specgram(x, Fs=fs)
axs[1].set_title("Noised")
axs[2].specgram(z, Fs=fs)
axs[2].set_title("Filtered")


axs[2].set(xlabel='Time s', ylabel='Hz')
# plt.savefig("data/speech/result.pdf")


import tikzplotlib
tikzplotlib.save("figures/tikz/spectrogram.tex")