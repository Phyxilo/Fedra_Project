//g++-5 eventSelector.cpp -w `root-config --cflags` -I$FEDRA_ROOT/include -L$FEDRA_ROOT/lib -lEIO -lDataConversion -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion `root-config --libs` -o eSel
//./eSel $1 $2 ... $n Enter event IDs as arguments

#include<stdio.h>
#include<EdbDataSet.h>
#include<vector>

using namespace std;

vector<string> allArgs;
vector<int> eventID;

int main(int argc, char *argv[])
{

	if (argc > 1)
	{
		allArgs.assign(argv + 1, argv + argc);

		for (int i = 0; i < argc-1; i++)
		{
			const char *eNum = allArgs[i].c_str();
			eventID.push_back(atoi(eNum));
			//cout << "Event ID: " << allArgs[i] << endl;
		}

		/*
		if (strstr(sentence.c_str(), word.c_str()))
		{
			cout << "OK" << endl;
		}
		*/
	}
	else{ cout << "Insert Event ID. Usage: ./eSel $eID" << endl; return 0; }

	EdbDataProc *dproc = new EdbDataProc("lnk.def");
	dproc->InitVolume( 100, "nseg>=3&&abs(t.eTX)<0.4&&abs(t.eTY)<0.4");
	EdbPVRec *pvr = dproc->PVR();

	int ntrk = pvr->Ntracks();

	TObjArray *eSelected = new TObjArray;

	for (int j = 0; j < eventID.size(); j++)
	{
		for (int i = 0; i < ntrk; i++)
		{
			EdbTrackP *t = pvr->GetTrack(i);
			//cout << t->MCEvt() << endl; //91150

			if (t->MCEvt() == eventID[j]) {eSelected->Add(t);}
		}

	}
	cout << "Track count of the event: " << eSelected->GetEntriesFast() << endl;

	dproc->MakeTracksTree(*eSelected, 0, 0, "selected.root");

	return 0;
}
