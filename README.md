# Morpha

Real-Time frame-by-frame morphing process implemented using spectral crossing.  
From: <https://www.dsprelated.com/freebooks/sasp/spectral_envelope_extraction.html>  

1. Short-Time Fourier Transform (STFT) of both the modulator and carrier signals;
2. Extraction of the spectral envelope of each time-frame of the signals;
3. Division of the spectrum of each carrier frame by its own spectral envelope in order to flattening it;
4. Multiplication of the flattened carrier spectral frame by the envelope of the corresponding modulator frame;
5. Inverse Short-Time Fourier Transform (ISTFT) of the resultant time-localized spectrum.

The code is based on the theory described in:  
[1] J. Smith. Spectral Audio Signal Processing.W3K Publishing, 2011  
[2] U. Zölzer. DAFX: Digital Audio Effects. Chichester, John Wiley & Sons, 2011  

Morphing factor controls how carrier, modulator or morphed signal is present.  
Using 0.0, only carrier is present, 0.5 morphed signal dominates and 1.0 only modulator is present.  
From 0.0 to 0.5 transition gradually shifts from carrier to morphed signal. From 0.5 to 1.0 it gradually transitions from morphed signal to modulator.  
cf value determines the cut-off frequency of the low-pass filter. A lower cf results in smoother spectral envelope shape.
