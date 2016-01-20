package Loading;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;

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

public class FFTbase implements Runnable {
	static double[]		realData;
	double[]			imagData;
	int					length	= 0;
	int					laenge;
	Complex[]			freq;
	Complex[]			oneMS;
	Complex[]			data;
	ArrayList<Double>	list;
	ArrayList<Double>	list2;
	int					bit;
	static int			n;
	float				samplesPerSec;

	static long			bytesPer10MS;

	FFTbase(int i, double[] realData) {
		this.laenge = i;
		this.realData = realData;
	}

	public Complex[] fft(Complex[] x) {
		int N = x.length;

		// base case
		if (N == 1)
			return new Complex[] {
				x[0]
			};

		// radix 2 Cooley-Tukey FFT
		int levels = 31 - Integer.numberOfLeadingZeros(N); // Equal to floor(log2(n))
		if (1 << levels != N)
			throw new IllegalArgumentException("Length is not a power of 2");

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

	public static void main(String[] args) throws UnsupportedAudioFileException, IOException {
		FFTbase fftb = new FFTbase(0, null);
		int r = 0;
		int n = 0;
		int l;
		realData = fftb.readFully(new File("10000.wav"));
		System.out.println("realdata: " + realData.length);
		double[] imagData = new double[realData.length];
		double[] range = new double[12];
		fftb.freq = new Complex[realData.length];

		// double[time in millisecs][frequenzen] 882bytes
		fftb.data = new Complex[(int) bytesPer10MS];
		// for (int i = 0; i < imagData.length; i++)
		// imagData[i] = 0;
		// for (int i = 0; i < realData.length; i++) {
		// data[i] = new Complex(realData[i], imagData[i]);
		// }

		// fft(data);
		Thread t1 = new Thread(new FFTbase(realData.length, realData));
		t1.start();

		l = fftb.retLength();
		double[] nfreq = new double[l / 12];
		while (r < l) {
			if (r + 12 < l) {
				for (int i = 0; i < 12; i++) {
					range[i] = fftb.freq[r].re;
					r++;
				}
			} else {
				break;
			}
			nfreq[n] = fftb.retMaxOfRange(range);
			n++;
		}

		System.out.println("hier kommt die ausgabe");
		// for (int j = 0; j < nfreq.length / 2; j++) {
		// System.out.println(j + ": " + nfreq[j]);
		// }
	}

	/**
	 * Berechnet ein Array welches 2^n groß ist und die Samples enthällt
	 * 
	 * @param file
	 * @return double[] samples
	 * @throws UnsupportedAudioFileException
	 * @throws IOException
	 */
	public double[] readFully(File file) throws UnsupportedAudioFileException, IOException {
		byte[] bytes = null;
		Double[] ret = null;
		Double[] samp = null;
		double[] samples = null;

		/*
		 * Einlesen der Ton-Datei
		 * Ausgeben von wichtigen Daten
		 */
		long time = -System.nanoTime();

		AudioInputStream audioInputStream = AudioSystem.getAudioInputStream(file);
		AudioFormat fmt = audioInputStream.getFormat();

		bit = fmt.getSampleSizeInBits();
		samplesPerSec = fmt.getSampleRate();
		bytesPer10MS = (long) (bit * samplesPerSec / 800);

		System.out.println(time + System.nanoTime() + "ns1");
		/*
		 * Audio wird zu Byte Array konvertiert
		 */
		long time2 = -System.nanoTime();
		try {
			if (fmt.getEncoding() != Encoding.PCM_SIGNED) {
				throw new UnsupportedAudioFileException();
			}
			// read the data fully
			bytes = new byte[audioInputStream.available()];
			audioInputStream.read(bytes);
		} finally {
			audioInputStream.close();
		}
		double max = Math.pow(2, bit - 1);

		ByteBuffer bb = ByteBuffer.wrap(bytes);
		bb.order(fmt.isBigEndian() ? ByteOrder.BIG_ENDIAN : ByteOrder.LITTLE_ENDIAN);

		samp = new Double[bytes.length * 8 / bit];
		System.out.println(time2 + System.nanoTime() + "ns2");
		/*
		 * convert sample-by-sample to a scale of
		 * -1.0 <= samples[i] < 1.0
		 * System.out.println(samp.length);
		 */
		long time3 = -System.nanoTime();
		for (int i = 0; i < samp.length; ++i) {
			switch (bit) {
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
		// for (double d : samp)
		// System.out.println(d);
		System.out.println(time3 + System.nanoTime() + "ns3");

		list = new ArrayList<Double>();

		long time4 = -System.nanoTime();
		// list = new ArrayList<Double>(Arrays.asList(samp));
		for (int i = 1; i < samp.length; i += 2) {
			list.add((Double) samp[i]);
		}
		System.out.println(time4 + System.nanoTime() + "ns4");
		System.out.println();
		// for (Double d : list)
		// System.out.println("List: " + d);

		/*
		 * berechnet die benötigte Größte des Input Array (muss 2^n sein, n€N)
		 */
		/*
		 * Erstellt ein neues Array mit der Größe 2^n und setzt die leeren Felder 0
		 */
		long time7 = -System.nanoTime();
		n = getPowerof2(list.size());
		samples = new double[n];
		for (int i = 0; i < n; i++) {
			if ((list.size()) > i)
				samples[i] = list.get(i);
			else
				samples[i] = 0;
		}
		System.out.println(time7 + System.nanoTime() + "ns7");
		// for (double d : samples)
		// System.out.println("Samples: " + d);
		return samples;
	}

	/**
	 * Ausgabe von RealData und ImagData
	 * 
	 * @param i
	 * @return String Complex
	 */
	public String toString(int i) {
		if (imagData[i] == 0)
			return Double.toString(realData[i]);
		if (realData[i] == 0)
			return imagData[i] + "i";
		if (imagData[i] < 0)
			return realData[i] + " - " + (-imagData[i]) + "i";
		return realData[i] + " + " + imagData[i] + "i";
	}

	/**
	 * Ausgabe von RealData
	 * 
	 * @param i
	 * @return realData
	 */
	public String toStringReal(int i) {
		return Double.toString(realData[i]);

	}

	/**
	 * Ausgabe von ImagData
	 * 
	 * @param i
	 * @return imagData
	 */
	public String toStringImag(int i) {
		return Double.toString(imagData[i]);

	}

	/**
	 * Ausgabe der Länge von freq
	 * 
	 * @return freq length
	 */
	public int retLength() {
		for (int j = 0; j < freq.length / 2; j++) {
			if (freq[j] == null)
				return j;
		}
		return freq.length;
	}

	/**
	 * Ausgabe des Maximalwertes von einem Array
	 * 
	 * @param range
	 * @return max
	 */
	public double retMaxOfRange(double[] range) {
		double max = 0;
		for (int i = 0; i < range.length; i++) {
			if (Math.abs(range[i]) > max)
				max = Math.abs(range[i]);
		}
		return max;
	}

	@Override
	public void run() {
		int r = 0;
		int l = laenge;
		int durchgang = 0;
		int bytes = getPowerof2((int) bytesPer10MS);
		data = new Complex[bytes];
		System.out.println(data.length);
		while (r < l) {
			long time4 = -System.nanoTime();
			if (r + (int) bytesPer10MS < l) {
				for (int i = 0; i < bytesPer10MS; i++, r++) {
					data[i] = new Complex(realData[r], 0);
				}
			} else
				break;
			for (int i = (int) bytesPer10MS; i < bytes; i++) {
				data[i] = new Complex(0, 0);
			}
			freq = fft(data);
			for (Complex c : freq)
				if (durchgang == 10)
					System.out.println(c.re);
			// System.out.println(time4 + System.nanoTime() + "ns run");
			try {
				durchgang++;
				// System.out.println("FERTIG!");
				Thread.sleep(5);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}

	/**
	 * Berechnung der nächst größeren 2er-Potenz von i
	 * 
	 * @param i
	 * @return nächst größeren 2er-Potenz von i
	 */
	public static int getPowerof2(int i) {
		int levels = 31 - Integer.numberOfLeadingZeros(i); // Equal to floor(log2(n))
		if (1 << levels != i)
			return Integer.highestOneBit(i * 2 - 1);
		return 1;
	}
}
