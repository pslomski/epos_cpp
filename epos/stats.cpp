// stats.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include "pch.h"
#include <iostream>
#include "epos.h"

using namespace std;
using namespace epos;

bool ProcessFile(const string& inFilename, const string& outFilename)
{
	epos::Channel channel;
	if (channel.Open(inFilename))
	{
		Histogram histogram(100, -1.0, 1.0);
		epos::Seg seg;
		size_t segCount = channel.GetSegCount();
		for (size_t segNr = 0; segNr < segCount; segNr++)
		{
			if (channel.ReadSeg(segNr, seg))
			{
				for (auto& val: seg.data)
					histogram.Add(val);
			}
			else
			{
				cout << "error reading segment" << endl;
				return false;
			}
		}
		histogram.SaveToFile(outFilename);
	}
	return true;
}

int main()
{
	string rootPath = "D:/_temp/epos/";
	ProcessFile(rootPath + "CeBr.bin", rootPath + "CeBr.hist");
	ProcessFile(rootPath + "RF.bin", rootPath + "RF.hist");
	ProcessFile(rootPath + "SiMP.bin", rootPath + "SiMP.hist");
}

