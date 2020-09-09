// epos2wav.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include "soundfile.h"
#include <stdlib.h>
#include <string>
#include <iostream>
#include "../epos/epos.h"

using namespace std;

bool ProcessFile(const string& inFilename, const string& outFilename)
{
	epos::Channel channel;
	if (channel.Open(inFilename))
	{
		SoundHeader header;
		header.setHighMono();
		SoundFileWrite soundfile(outFilename.c_str(), header);
		epos::Seg seg;
		size_t segCount = channel.GetSegCount();
		for (size_t segNr = 0; segNr < segCount; segNr++)
		{
			if (channel.ReadSeg(segNr, seg))
			{
				for (size_t i = 0; i < seg.data.size(); i++)
				{
					soundfile.writeSampleDouble(seg.data[i]);
				}
			}
			else
			{
				cout << "error reading segment" << endl;
				return false;
			}
		}
		// file will be closed automatically, but can be done manually:
		soundfile.close();
	}
	return true;
}

int main(int argc, char** argv)
{
	string rootPath = "D:/_temp/epos/";
	ProcessFile(rootPath + "CeBr.bin", rootPath + "CeBr.wav");
	ProcessFile(rootPath + "RF.bin", rootPath + "RF.wav");
	ProcessFile(rootPath + "SiMP.bin", rootPath + "SiMP.wav");


	return 0;
}