#ifndef LIB_EPOS_HPP_INCLUDED
#define LIB_EPOS_HPP_INCLUDED

#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace epos
{
	class EEposError : public std::runtime_error {

	};

	typedef double Float;
	typedef short int Sample; // Raw data type
	typedef std::vector<int> VectorInt;
	typedef std::vector<Float> VectorFloat;
	typedef VectorFloat::iterator VectorFloatIter;
	typedef std::vector<Sample> SegBuf; // Segment buffer data
	typedef std::string String;
	typedef const std::string& StringRef;

	constexpr Float time_resolution = 0.5; // Duration of sample in nanosec [ns]
	constexpr size_t NbrSAMP = 1600; // Number of samples in segment
	constexpr size_t NbrSEG = 1000; // Number of segments in single data file

	template <typename T> int sgn(T val)
	{
		return (T(0) < val) - (val < T(0));
	}

	template <typename T> bool isInRange(T val, T low, T high)
	{
		return (low <= val && val <= high);
	}

	template <typename T> bool are_equal(T expected, T actual, T precision = 1e-6)
	{
		return abs(expected - actual) <= precision;
	}

	// Arithmetic mean of [begin, end) range
	Float arithmetic_mean(const VectorFloatIter begin, const VectorFloatIter end);

	// Standard deviation of [begin, end) range
	Float variance(const VectorFloatIter begin, const VectorFloatIter end);

	void avevar(VectorFloatIter begin, VectorFloatIter end, Float& ave, Float& var);

	Float nachylenie5(VectorFloat& data, int i0);

	struct SeriesStatistics {
		Float arithmetic_mean;
		Float variance;
		Float standard_deviation;
		void compute_from_range(VectorFloatIter begin, VectorFloatIter end);
	};

	class LookupTable {
	private:
		long long tableSize;
		Float A;
		Float B;
		Float indexScale;
		VectorFloat ftable;

		union InToRFloat {
			double f;
			long long i;
		};

		void initftable(Float f(Float))
		{
			ftable.resize(tableSize);
			for (int i = 0; i < tableSize; ++i)
				ftable[i] = f((Float)i / indexScale + A);
		}

	public:
		LookupTable(long long _tableSize, Float _A, Float _B, Float f(Float)): tableSize(_tableSize), A(_A), B(_B)
		{
			indexScale = (Float)tableSize / (B - A);
			initftable(f);
		}

		Float flut(Float x)
		{
			/*constexpr Float ftoiBias = 12582912.0;
			InToRFloat fTmp;
			fTmp.f = x * indexScale + (ftoiBias - A * indexScale);
			long long i = fTmp.i & (tableSize - 1);
			return ftable[i];*/
			//return ftable[x * (indexScale + A)];
			return 0.0;
		}
	};

	// Class representing data segment of single channel
	class Channel {
	public:
		class Stats {
		public:
			Float DC; // Background level of the channel
			Float DCstddev; // Standard deviation of the background level
							// W dwa odchylenia standardowe wpada 0,95449989 wartosci.
			Float min; // Minimum sample value 
			Float max; // Maximum sample value 
			int idxmin; // Index of minimum sample value
			int idxmax; // Index of maximum sample value
			void compute_statistics(const VectorFloat& data);
		};

		class Result {
		public:
			Float t0;
			VectorFloat a;
		};

	public:
		static constexpr double dt = 0.5; // time resolution

		// sampleIndex from sampleTime
		// for sampleIndex = 0, T = 0.5
		// for sampleIndex = NbrSAMP -1, T = 800.0
		static double getT(size_t sampleIndex) {
			return (sampleIndex + 1) * dt;
		}

		// sampleTime from sampleIndex
		// dla sampleTime = 0.5, Idx = 0
		// dla sampleTime = 800.0, Idx = NbrSAMP -1
		static size_t getIdx(double sampleTime) {
			constexpr double dt_inv = 1.0 / Channel::dt;
			return static_cast<int>(sampleTime * dt_inv - 1);
		}

		static VectorFloat time; // time domain
		VectorFloat data;
		Stats stats;
		Result res;

		void init(size_t dataSize);
		void add(Float val) {
			data.push_back(val);
		}
		Stats& computeStats();
		double getMean(size_t SampleCount) const;
		//void addConst(Float c);
		void saveToFile(StringRef filePath) const;
	};

	// Klasa Channel reprezentuje kana³ danych w formacie binarnym
	class BinChannel {
	private:
		// rozmiar segmentu w bajtach
		static constexpr size_t bufByteCount = epos::NbrSAMP * sizeof(Sample);
		static std::streampos getSegStartPos(size_t segNbr) {
			return segNbr * BinChannel::bufByteCount;
		}

		std::streamsize fFileSize;
		std::ifstream fFile; // aktualnie przetwarzany plik
	public:
		BinChannel();
		~BinChannel();
		// Opens epos data file
		bool open(StringRef filePath);
		// Closes epos data file
		void close();
		// Read segment of number segNbr
		// Segments numbers range from 0 to getSegCount - 1
		bool read(size_t segNr, Channel& seg);
		size_t getSegCount() const {
			return static_cast<size_t>(fFileSize / bufByteCount);
		}
	};

	class Segment {
	public:
		size_t fileNr; // numer pliku
		size_t segNr; // Segment number from 0 to epos::NbrSEG - 1
		bool clipped = false;
		Channel ch1;
		Channel ch2;
		Channel ch3;
		void saveToFile(StringRef filePath) const;
	};

	// Interface for dataset manipulation
	class DataSet {
	public:
		virtual bool read(Segment& data) = 0;
	};

	//klasa odczytujaca kanaly danych w formacie binarnym
	class DataSetBin : public DataSet {
	private:
		size_t fCurrSeg;
		size_t fSegMax;
		BinChannel fCh1;
		BinChannel fCh2;
		BinChannel fCh3;
		bool isValid() const;
	public:
		bool open(StringRef rootPath, size_t segMin = 0, size_t segMax = std::numeric_limits<short int>::max());
		bool read(Segment& data);
		bool getData(size_t segNr, Segment& data);
	};

	//klasa odczytujaca dane z plikow tekstowych
	class DataSetTxt : public DataSet {
	private:
		String fRootPath;
		std::ifstream fDataFile;
		Float fScalingFactor;
		int fCurrent;
		int fStop;
		int fSegNr;

		bool openFile(int fileNr);
	public:
		// Samples are multiplied by scalingFactor
		DataSetTxt(Float scalingFactor = 1.0);
		bool open(StringRef rootPath, int fileMin, int fileMax);
		bool read(Segment& data);
	};

	// Class for reading raw short int samples
	class DataSetTxtRaw : public DataSetTxt {
	public:
		DataSetTxtRaw() :
			DataSetTxt(1.0) {}
	};

	// Class for reading samples normalized to (-1..1)
	class DataSetTxtNorm : public DataSetTxt {
	public:
		DataSetTxtNorm() :
			DataSetTxt(1.0 / std::numeric_limits<short int>::max()) {}
	};

	class Filter {
	protected:
		double fAlpha;
		double fXprev;
		double fYprev;
	public:
		virtual void Init(double RC, double dt, double startVal) = 0;
		virtual double Get(double x) = 0;
		virtual void FilterBuf(double RC, double dt, const VectorFloat& in, VectorFloat& out) = 0;
	};

	class Highpass : public Filter {
	private:
	public:
		void Init(double RC, double dt, double startVal);
		double Get(double x);
		void FilterBuf(double RC, double dt, const VectorFloat& in, VectorFloat& out);
	};

	class Lowpass : public Filter {
	private:
	public:
		void Init(double RC, double dt, double startVal);
		double Get(double x);
		void FilterBuf(double RC, double dt, const VectorFloat& in, VectorFloat& out);
	};

	typedef size_t PeakIndex;

	class Peak {
	public:
		Float t0;
		Float y0;
		PeakIndex t0idx;
		Float A;
	};

	typedef std::vector<Peak> Peaks;

	class PeakDetector {
	public:
		static Peaks get_all(PeakDetector& peak_detector);

		virtual ~PeakDetector() {}
		virtual bool next(Peak &peak) = 0;
	};

} // namespace epos

#endif // LIB_EPOS_HPP_INCLUDED