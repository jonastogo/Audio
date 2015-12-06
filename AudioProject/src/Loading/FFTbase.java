package Loading;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioFormat.Encoding;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.UnsupportedAudioFileException;

//Java Implementierung der Fast Fourier Transformation
//Braucht Complex.java aus der "Einfuehrung in der Rechnerbedienung"
// http://www.theorie.physik.uni-goettingen.de/~honecker/rb07/Complex.java

//*** Die eigentliche FFT
//f[]: Eingabe: zu transformierende Daten
//   Ausgabe: Ergebnis der Transformation
//sign=-1: Hintransformation; sign=1: Ruecktransformation

public class FFTbase {
	static double[]		realData;
	static double[]		imagData;
	static int			length	= 0;
	static Complex[]	freq;

	public static Complex[] fft(Complex[] x) {
		int N = x.length;

		// base case
		if (N == 1)
			return new Complex[] {
				x[0]
			};

		// radix 2 Cooley-Tukey FFT
		if (N % 2 != 0) {
			throw new RuntimeException("N is not a power of 2");
		}

		// fft of even terms
		Complex[] even = new Complex[N / 2];
		for (int k = 0; k < N / 2; k++) {
			even[k] = x[2 * k];
		}
		Complex[] q = fft(even);

		// fft of odd terms
		Complex[] odd = even; // reuse the array
		for (int k = 0; k < N / 2; k++) {
			odd[k] = x[2 * k + 1];
		}
		Complex[] r = fft(odd);

		// combine
		Complex[] y = new Complex[N];
		for (int k = 0; k < N / 2; k++) {
			double kth = -2 * k * Math.PI / N;
			Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
			y[k] = q[k].plus(wk.times(r[k]));
			y[k + N / 2] = q[k].minus(wk.times(r[k]));
		}
		freq = y;
		return y;
	}

	/*
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function.
	 */
	public static void transform(double[] real, double[] imag) {
		if (real.length != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");

		int n = real.length;
		if (n == 0)
			return;
		else if ((n & (n - 1)) == 0) // Is power of 2
			transformRadix2(real, imag);
		else
			// More complicated algorithm for arbitrary sizes
			transformBluestein(real, imag);
	}

	/*
	 * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function. This transform does not perform scaling, so the inverse is not a true inverse.
	 */
	public static void inverseTransform(double[] real, double[] imag) {
		transform(imag, real);
	}

	/*
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
	 */
	public static void transformRadix2(double[] real, double[] imag) {
		// Initialization
		if (real.length != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		int n = real.length;
		int levels = 31 - Integer.numberOfLeadingZeros(n); // Equal to floor(log2(n))
		if (1 << levels != n)
			throw new IllegalArgumentException("Length is not a power of 2");
		double[] cosTable = new double[n / 2];
		double[] sinTable = new double[n / 2];
		for (int i = 0; i < n / 2; i++) {
			cosTable[i] = Math.cos(2 * Math.PI * i / n);
			sinTable[i] = Math.sin(2 * Math.PI * i / n);
		}

		// Bit-reversed addressing permutation
		for (int i = 0; i < n; i++) {
			int j = Integer.reverse(i) >>> (32 - levels);
			if (j > i) {
				double temp = real[i];
				real[i] = real[j];
				real[j] = temp;
				temp = imag[i];
				imag[i] = imag[j];
				imag[j] = temp;
			}
		}

		// Cooley-Tukey decimation-in-time radix-2 FFT
		for (int size = 2; size <= n; size *= 2) {
			int halfsize = size / 2;
			int tablestep = n / size;
			for (int i = 0; i < n; i += size) {
				for (int j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
					double tpre = real[j + halfsize] * cosTable[k] + imag[j + halfsize] * sinTable[k];
					double tpim = -real[j + halfsize] * sinTable[k] + imag[j + halfsize] * cosTable[k];
					real[j + halfsize] = real[j] - tpre;
					imag[j + halfsize] = imag[j] - tpim;
					real[j] += tpre;
					imag[j] += tpim;
				}
			}
			if (size == n) // Prevent overflow in 'size *= 2'
				break;
		}
		realData = real;
		imagData = imag;
	}

	/*
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
	 * Uses Bluestein's chirp z-transform algorithm.
	 */
	public static void transformBluestein(double[] real, double[] imag) {
		// Find a power-of-2 convolution length m such that m >= n * 2 + 1
		if (real.length != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		int n = real.length;
		if (n >= 0x20000000)
			throw new IllegalArgumentException("Array too large");
		int m = Integer.highestOneBit(n * 2 + 1) << 1;

		// Trignometric tables
		double[] cosTable = new double[n];
		double[] sinTable = new double[n];
		for (int i = 0; i < n; i++) {
			int j = (int) ((long) i * i % (n * 2)); // This is more accurate than j = i * i
			cosTable[i] = Math.cos(Math.PI * j / n);
			sinTable[i] = Math.sin(Math.PI * j / n);
		}

		// Temporary vectors and preprocessing
		double[] areal = new double[m];
		double[] aimag = new double[m];
		for (int i = 0; i < n; i++) {
			areal[i] = real[i] * cosTable[i] + imag[i] * sinTable[i];
			aimag[i] = -real[i] * sinTable[i] + imag[i] * cosTable[i];
		}
		double[] breal = new double[m];
		double[] bimag = new double[m];
		breal[0] = cosTable[0];
		bimag[0] = sinTable[0];
		for (int i = 1; i < n; i++) {
			breal[i] = breal[m - i] = cosTable[i];
			bimag[i] = bimag[m - i] = sinTable[i];
		}

		// Convolution
		double[] creal = new double[m];
		double[] cimag = new double[m];
		convolve(areal, aimag, breal, bimag, creal, cimag);

		// Postprocessing
		for (int i = 0; i < n; i++) {
			real[i] = creal[i] * cosTable[i] + cimag[i] * sinTable[i];
			imag[i] = -creal[i] * sinTable[i] + cimag[i] * cosTable[i];
		}
	}

	/*
	 * Computes the circular convolution of the given real vectors. Each vector's length must be the same.
	 */
	public static void convolve(double[] x, double[] y, double[] out) {
		if (x.length != y.length || x.length != out.length)
			throw new IllegalArgumentException("Mismatched lengths");
		int n = x.length;
		convolve(x, new double[n], y, new double[n], out, new double[n]);
	}

	/*
	 * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
	 */
	public static void convolve(double[] xreal, double[] ximag, double[] yreal, double[] yimag, double[] outreal, double[] outimag) {
		if (xreal.length != ximag.length || xreal.length != yreal.length || yreal.length != yimag.length || xreal.length != outreal.length || outreal.length != outimag.length)
			throw new IllegalArgumentException("Mismatched lengths");

		int n = xreal.length;
		xreal = xreal.clone();
		ximag = ximag.clone();
		yreal = yreal.clone();
		yimag = yimag.clone();

		transform(xreal, ximag);
		transform(yreal, yimag);
		for (int i = 0; i < n; i++) {
			double temp = xreal[i] * yreal[i] - ximag[i] * yimag[i];
			ximag[i] = ximag[i] * yreal[i] + xreal[i] * yimag[i];
			xreal[i] = temp;
		}
		inverseTransform(xreal, ximag);
		for (int i = 0; i < n; i++) { // Scaling (because this FFT implementation omits it)
			outreal[i] = xreal[i] / n;
			outimag[i] = ximag[i] / n;
		}
	}

	public static void main(String[] args) throws UnsupportedAudioFileException, IOException {
		int r = 0;
		int n = 0;
		int l;
		double[] realData = readFully(new File("Windows Logon.wav"));
		double[] imagData = new double[realData.length];
		double[] range = new double[12];

		// double[time in millisecs][frequenzen] 882bytes
		Complex[] data = new Complex[realData.length];
		freq = new Complex[realData.length];
		for (int i = 0; i < imagData.length; i++)
			imagData[i] = 0;
		for (int i = 0; i < realData.length; i++) {
			data[i] = new Complex(realData[i], imagData[i]);
		}

		fft(data);

		l = retLength();
		System.out.println("laenge: " + l);
		double[] nfreq = new double[l / 12];
		while (r < l) {
			if (r + 12 < l) {
				for (int i = 0; i < 12; i++) {
					range[i] = freq[r].re;
					r++;
				}
			} else {
				System.out.println("hallo?");
				break;
			}
			nfreq[n] = retMaxOfRange(range);
			n++;
		}

		System.out.println("hier kommt die ausgabe");
		for (int j = 0; j < nfreq.length / 2; j++) {
			System.out.println(nfreq[j]);
		}
	}

	public static double[] readFully(File file) throws UnsupportedAudioFileException, IOException {
		byte[] bytes = null;
		double[] ret = null;
		double[] samp = null;
		double[] samples = null;
		AudioInputStream audioInputStream = AudioSystem.getAudioInputStream(file);
		AudioFormat fmt = audioInputStream.getFormat();
		try {
			if (fmt.getEncoding() != Encoding.PCM_SIGNED) {
				throw new UnsupportedAudioFileException();
			}

			// read the data fully
			bytes = new byte[audioInputStream.available()];
			System.out.println(audioInputStream.available());
			audioInputStream.read(bytes);
		} finally {
			audioInputStream.close();
		}

		int bits = fmt.getSampleSizeInBits();
		System.out.println("bits: " + bits);
		double max = Math.pow(2, bits - 1);

		ByteBuffer bb = ByteBuffer.wrap(bytes);
		bb.order(fmt.isBigEndian() ? ByteOrder.BIG_ENDIAN : ByteOrder.LITTLE_ENDIAN);

		samp = new double[bytes.length * 8 / bits];
		// convert sample-by-sample to a scale of
		// -1.0 <= samples[i] < 1.0
		// System.out.println(samp.length);
		for (int i = 0; i < samp.length; ++i) {
			switch (bits) {
				case 8:
					samp[i] = (bb.get() / max);
					// System.out.println(bb.get() / max);
					break;
				case 16:
					samp[i] = (bb.getShort() / max);
					// System.out.println(bb.getShort() / max);
					break;
				case 32:
					samp[i] = (bb.getInt() / max);
					// System.out.println(bb.getInt() / max);
					break;
				case 64:
					samp[i] = (bb.getLong() / max);
					// System.out.println(bb.getLong() / max);
					break;
				default:
					throw new UnsupportedAudioFileException();
			}
		}

		ret = new double[(samp.length / 2)];
		length = ret.length;
		System.out.println("retlänge: " + ret.length + "  samplänge: " + samp.length + "  byteslänge: " + bytes.length);
		for (int i = 1, j = 0; i < samp.length; i += 2, j++)
			ret[j] = samp[i];

		int n = 1;
		while (true) {
			long length = (long) Math.pow(2, n);
			if (length >= ret.length)
				break;
			n += 1;
		}
		long samplesLength = (long) Math.pow(2, n);
		samples = new double[(int) samplesLength];
		System.out.println(samples.length);
		for (int i = 0; i < samplesLength; i++) {
			if ((ret.length) > i)
				samples[i] = ret[i];
			else
				samples[i] = 0;
		}
		return samples;
	}

	public static String toString(int i) {
		if (imagData[i] == 0)
			return Double.toString(realData[i]);
		if (realData[i] == 0)
			return imagData[i] + "i";
		if (imagData[i] < 0)
			return realData[i] + " - " + (-imagData[i]) + "i";
		return realData[i] + " + " + imagData[i] + "i";
	}

	public static String toStringReal(int i) {
		return Double.toString(realData[i]);

	}

	public static String toStringImag(int i) {
		return Double.toString(imagData[i]);

	}

	public static int retLength() {
		for (int j = 0; j < freq.length / 2; j++) {
			if (freq[j] == null)
				return j;
		}
		return freq.length;
	}

	public static double retMaxOfRange(double[] range) {
		double max = 0;
		for (int i = 0; i < range.length; i++) {
			if (Math.abs(range[i]) > max)
				max = Math.abs(range[i]);
		}
		return max;
	}
}
