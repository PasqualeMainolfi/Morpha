#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <array>
#include <numeric>
#include <vector>
#include <sndfile.h>
#include <fftw3.h>

using cframe = std::vector<std::array<double, 2>>;
using dframe = std::vector<double>;

/*
1. Short-Time Fourier Transform (STFT) of both the modulator and carrier signals;
2. Extraction of the spectral envelope of each time-frame of the signals;
3. Division of the spectrum of each carrier frame by its own spectral envelope in order to flattening it;
4. Multiplication of the flattened carrier spectral frame by the envelope of the corresponding modulator frame;
5. Inverse Short-Time Fourier Transform (ISTFT) of the resultant time-localized spectrum.

The code is based on the theory described in:
	[1] J. Smith. Spectral Audio Signal Processing.W3K Publishing, 2011.
	[2] U. Zölzer. DAFX: Digital Audio Effects. Chichester, John Wiley & Sons, 2011.
*/


class Morpha
{
private:
	size_t buffer_length;
	size_t half_buffer;
	fftw_plan fft_planner;
	fftw_plan ifft_planner;
	double* fftin;
	fftw_complex* fftout;
	fftw_complex* ifftin;
	double* ifftout;

public:
	Morpha(size_t buffer_length)
	:
	buffer_length(buffer_length),
	half_buffer(buffer_length / 2 + 1)
	{
		this->fftin = (double*) fftw_malloc(sizeof(double) * this->buffer_length);
		this->fftout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (this->half_buffer));
		this->fft_planner = fftw_plan_dft_r2c_1d(this->buffer_length, this->fftin, this->fftout, FFTW_MEASURE);

		this->ifftin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->half_buffer);
		this->ifftout = (double*) fftw_malloc(sizeof(double) * this->buffer_length);
		this->ifft_planner = fftw_plan_dft_c2r_1d(this->buffer_length, this->ifftin, this->ifftout, FFTW_MEASURE);

	}

	~Morpha()
	{
		fftw_destroy_plan(this->fft_planner);
		fftw_destroy_plan(this->ifft_planner);
		free(this->fftin);
		free(this->fftout);
		free(this->ifftin);
		free(this->ifftout);
	}

	void fft(dframe& source) {
		std::memcpy(this->fftin, source.data(), sizeof(double) * this->buffer_length);
		fftw_execute(this->fft_planner);
	}

	dframe ifft(cframe& data, bool norm) {
		std::memcpy(this->ifftin, data.data(), sizeof(fftw_complex) * this->half_buffer);
		fftw_execute(this->ifft_planner);
		dframe y(this->buffer_length, 0.0);
		for (size_t i = 0; i < this->buffer_length; ++i) {
			y[i] = norm == true ? this->ifftout[i] / (double) this->buffer_length : this->ifftout[i];
		}
		return y;
	}

	dframe get_spectral_envelope(dframe& data, size_t cf, bool scale_factor) {

		this->fft(data);

		dframe source_mag(this->half_buffer);
		cframe source_log(this->half_buffer);
		for (size_t i = 0; i < this->half_buffer; ++i) {
			double re = this->fftout[i][0];
			double im = this->fftout[i][1];
			double mag = std::sqrt(re * re + im * im);
			source_mag[i] = mag;
			source_log[i] = { std::log(mag + 1e-12), 0.0 };
		}

		dframe ifft_log = this->ifft(source_log, true);
		for (size_t i = 0; i < this->half_buffer; ++i) {
			double w = i == cf ? 0.5 : (i < cf ? 1.0 : 0.0);
			ifft_log[i] *= w;
		}

		this->fft(ifft_log);

		dframe real_cepstrum(this->half_buffer, 0.0);
		for (size_t i = 0; i < this->half_buffer; ++i) {
			real_cepstrum[i] = this->fftout[i][0];
		}

		double mean = std::accumulate(real_cepstrum.begin(), real_cepstrum.end(), 0.0) / this->buffer_length;
		for (size_t i = 0; i < this->half_buffer; ++i) {
			real_cepstrum[i] = std::exp(real_cepstrum[i] - mean);
		}

		double max_mag = *std::max_element(source_mag.begin(), source_mag.end());
		double max_cep = *std::max_element(real_cepstrum.begin(), real_cepstrum.end());
		double scale = max_mag / max_cep;

		if (scale_factor) {
			for (size_t i = 0; i < this->buffer_length; ++i) {
				real_cepstrum[i] *= scale;
			}
		}

		return real_cepstrum;
	}

	dframe morphing_process(dframe& carrier, dframe& modulator, size_t cf, double morphing_factor) {

		dframe source_cepstrum = this->get_spectral_envelope(carrier, cf, true);
		dframe target_cepstrum = this->get_spectral_envelope(modulator, cf, true);

		this->fft(modulator);

		cframe target_flatten_spectrum(this->half_buffer);
		for (size_t i = 0; i < this->half_buffer; ++i) {
			double re = this->fftout[i][0] / target_cepstrum[i];
			double im = this->fftout[i][1] / target_cepstrum[i];
			target_flatten_spectrum[i] = { re, im };
		}

		cframe morphed(this->half_buffer);
		for (size_t i = 0; i < this->half_buffer; ++i) {
			double re = source_cepstrum[i] * target_flatten_spectrum[i][0];
			double im = source_cepstrum[i] * target_flatten_spectrum[i][1];
			morphed[i] = { re, im };
		}

		this->fft(carrier);
		cframe carrier_fft(this->half_buffer, { 0.0, 0.0 });
		for (size_t i = 0; i < this->half_buffer; ++i) {
			carrier_fft[i][0] = this->fftout[i][0];
			carrier_fft[i][1] = this->fftout[i][1];
		}

		this->fft(modulator);
		cframe modulator_fft(this->half_buffer, { 0.0, 0.0 });
		for (size_t i = 0; i < this->half_buffer; ++i) {
			modulator_fft[i][0] = this->fftout[i][0];
			modulator_fft[i][1] = this->fftout[i][1];
		}

		double sx, cx, dx;

		if (morphing_factor <= 0.5) {
			sx = 1.0 - morphing_factor;
			cx = 2.0 * morphing_factor;
			dx = 0.0;
		} else {
			sx = 0.0;
			cx = 2.0 * (1.0 - morphing_factor);
			dx = 2.0 * (morphing_factor - 0.5);
		}

		cframe result(this->half_buffer, { 0.0, 0.0 });
		for (size_t i = 0; i < this->half_buffer; ++i) {
			double re = sx * carrier_fft[i][0] + cx * morphed[i][0] + dx * modulator_fft[i][0];
			double im = sx * carrier_fft[i][1] + cx * morphed[i][1] + dx * modulator_fft[i][1];
			result[i][0] = re / 2.0;
			result[i][1] = im / 2.0;
		}

		dframe m = this->ifft(result, false);
		double max_morphed = 0.0;
		for (size_t i = 0; i < this->buffer_length; ++i) {
			double current_value = std::abs(m[i]);
			max_morphed = current_value > max_morphed ? current_value : max_morphed;
		}

		for (size_t i = 0; i < this->buffer_length; ++i) {
			m[i] /= max_morphed;
		}

		return m;

	}
};


int main() {

	size_t buffer_length = 88200;

	SF_INFO info_source;
	info_source.format = 0;
	SNDFILE* source = sf_open("/Users/pm/AcaHub/AudioSamples/carrier2244100.wav", SFM_READ, &info_source);
	sf_count_t source_nsamples = info_source.frames * info_source.channels;
	dframe source_vector(source_nsamples, 0.0);
	sf_read_double(source, source_vector.data(), source_nsamples);

	SF_INFO info_target;
	info_target.format = 0;
	SNDFILE* target = sf_open("/Users/pm/AcaHub/AudioSamples/modulator2244100.wav", SFM_READ, &info_target);
	sf_count_t target_nsamples = info_target.frames * info_target.channels;
	dframe target_vector(target_nsamples, 0.0);
	sf_read_double(target, target_vector.data(), target_nsamples);

	dframe source_frame(buffer_length, 0.0);
	dframe target_frame(buffer_length, 0.0);

	for (size_t i = 0; i < buffer_length; ++i) {
		source_frame[i] = source_vector[i];
		target_frame[i] = target_vector[i];
	}

	sf_close(source);
	sf_close(target);

	Morpha m(buffer_length);
	dframe morphed = m.morphing_process(source_frame, target_frame, 22050, 0.5);

	SF_INFO info_morphed;
	info_morphed.samplerate = 44100;
	info_morphed.channels = 1;
	info_morphed.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;

	SNDFILE* audio_morphed = sf_open("test_morphed.wav", SFM_WRITE, &info_morphed);
	sf_write_double(audio_morphed, morphed.data(), morphed.size());

	sf_close(audio_morphed);

    return 0;
}
