//g++-5 dTrCut.cpp -w `root-config --cflags` -I$FEDRA_ROOT/include -L$FEDRA_ROOT/lib -lEIO -lDataConversion -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion `root-config --libs` -o dTrCut
//./dTrCut

#include<stdio.h>
#include<EdbDataSet.h>
#include<list>

using namespace std;

EdbTrackP *priTr = new EdbTrackP(8);
EdbTrackP *secTr = new EdbTrackP(8);

EdbTrackP *ObjArrCompare (TObjArray *priArr, TObjArray *secArr, int index);

int main(int argc, char *argv[]){

	// Declear the EdbDataProc object with the definition file "lnk.def"
	EdbDataProc *dproc = new EdbDataProc("lnk.def");

	// Read track data (data type=100 means read only the reconstructed tracks from linked_tracks.root)
	dproc->InitVolume( 100, "nseg>=3&&abs(t.eTX)<0.4&&abs(t.eTY)<0.4");

	// Get EdbPVRec object
	EdbPVRec *pvr = dproc->PVR();

	// Loop over the tracks
	int ntrk = pvr->Ntracks();

	TObjArray *slpSel = new TObjArray;
	TObjArray *effSel = new TObjArray;

	Float_t trTX;
	Float_t trTY;

	float slpCut = 0.4;

	int refIDEvent = 0;
	list<int> eIDList;

	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);

		slpSel->Add(t);
	}
	cout << "Selected tracks count after small angle cut: " << slpSel->GetEntriesFast() << endl;

	EdbTrackP *refTrack = new EdbTrackP(8);
	EdbTrackP *secTrack = new EdbTrackP(8);

	Float_t refTrTX; Float_t refTrTY;
	Float_t secTrTX; Float_t secTrTY;

	Float_t dTrTX; Float_t dTrTY;

	for(int i = 0; i < slpSel->GetEntriesFast(); i++)
	{
		refTrack = (EdbTrackP*)(slpSel->At(i));
		refTrTX = refTrack->TX(); refTrTY = refTrack->TY();
		refIDEvent = refTrack->MCEvt();

		if (effSel->FindObject(refTrack) == 0 && find(eIDList.begin(), eIDList.end(), refIDEvent) != eIDList.end())
		{
			for(int j = 0; j < slpSel->GetEntriesFast(); j++)
			{
				secTrack = (EdbTrackP*)(slpSel->At(j));
				secTrTX = secTrack->TX(); secTrTY = secTrack->TY();

				if (refTrack != secTrack)
				{
					dTrTX = refTrTX - secTrTX;
					dTrTY = refTrTY - secTrTY;

					/*
					if (refIDEvent == secTrack->MCEvt() && (dTrTX > 0.02 || dTrTX < -0.02) && (dTrTY > 0.02 || dTrTY < -0.02) && effSel->FindObject(secTrack) == 0)
					{
						effSel->Add(secTrack);
					}
					*/
					if (refIDEvent == secTrack->MCEvt() && !(dTrTX < 0.02 && dTrTX > -0.02 && dTrTY < 0.02 && dTrTY > -0.02) && effSel->FindObject(secTrack) == 0)
					{
						effSel->Add(secTrack);
					}
				}
			}
		}
		eIDList.push_back(refIDEvent);
	}
	cout << "Selected tracks count: " << effSel->GetEntriesFast() << endl;

	dproc->MakeTracksTree(*slpSel, 0, 0, "Out/primaryCut.root");
	dproc->MakeTracksTree(*effSel, 0, 0, "Out/finalCut.root");

	return 0;
}
