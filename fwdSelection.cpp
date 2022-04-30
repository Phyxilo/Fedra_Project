//g++-5 fwdSelection.cpp -w `root-config --cflags` -I$FEDRA_ROOT/include -L$FEDRA_ROOT/lib -lEIO -lDataConversion -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion `root-config --libs` -o fwdSel
//./fwdSel lnk.def

#include<stdio.h>
#include<EdbDataSet.h>
#include<list>

using namespace std;

EdbTrackP *cmpTr = new EdbTrackP(8);
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
	TObjArray *dscSel = new TObjArray;

	Float_t trTX;
	Float_t trTY;

	float slpCut = 0.4;

	int refIDEvent = 0;
	list<int> eIDList;

	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);
		trTX = t->TX(); trTY = t->TY();
		if (trTX <= slpCut && trTX >= -slpCut && trTY <= slpCut && trTY >= -slpCut)
		{
			slpSel->Add(t);
		}
	}
	cout << "Selected tracks count after small angle cut: " << slpSel->GetEntriesFast() << endl;

	EdbTrackP *refTrack = new EdbTrackP(8);
	EdbTrackP *secTrack = new EdbTrackP(8);

	Float_t refTrX; Float_t refTrY;
	Float_t secTrX; Float_t secTrY;
	Float_t refTrTX; Float_t refTrTY;
	Float_t secTrTX; Float_t secTrTY;

	Float_t dTrX; Float_t dTrY;
	Float_t dTrTX; Float_t dTrTY;

	Float_t trDist;

	for(int i = 0; i < slpSel->GetEntriesFast(); i++)
	{
		refTrack = (EdbTrackP*)(slpSel->At(i));
		refTrX = refTrack->X(); refTrY = refTrack->Y();
		refTrTX = refTrack->TX(); refTrTY = refTrack->TY();
		refIDEvent = refTrack->MCEvt();

		if (dscSel->FindObject(refTrack) == 0 && find(eIDList.begin(), eIDList.end(), refIDEvent) != eIDList.end())
		{
			for(int j = 0; j < slpSel->GetEntriesFast(); j++)
			{
				secTrack = (EdbTrackP*)(slpSel->At(j));
				secTrX = secTrack->X(); secTrY = secTrack->Y();
				secTrTX = secTrack->TX(); secTrTY = secTrack->TY();

				if (refTrack != secTrack)
				{
					dTrX = refTrX - secTrX;
					dTrY = refTrY - secTrY;

					dTrTX = refTrTX - secTrTX;
					dTrTY = refTrTY - secTrTY;

					trDist = sqrt(dTrX*dTrX + dTrY*dTrY);

					if (refIDEvent == secTrack->MCEvt() && dTrTX < 0.02 && dTrTX > -0.02 && dTrTY < 0.02 && dTrTY > -0.02 && dscSel->FindObject(secTrack) == 0)
					{
						dscSel->Add(secTrack);
					}
				}
			}
		}
		eIDList.push_back(refIDEvent);
	}
	cout << "Discarded tracks count: " << dscSel->GetEntriesFast() << endl;

	for (int i = 0; i < slpSel->GetEntriesFast(); i++)
	{
		cmpTr = ObjArrCompare(slpSel, dscSel, i);
		if (cmpTr != 0) {effSel->Add(cmpTr);}
	}
	cout << effSel->GetEntriesFast() << endl;

	dproc->MakeTracksTree(*slpSel, 0, 0, "Out/primaryCut.root");
	dproc->MakeTracksTree(*dscSel, 0, 0, "Out/discarded.root");
	dproc->MakeTracksTree(*effSel, 0, 0, "Out/finalCut.root");

	return 0;
}

EdbTrackP *ObjArrCompare (TObjArray *priArr, TObjArray *secArr, int index)
{
	priTr = (EdbTrackP*)(priArr->At(index));

	for (int j = 0; j < secArr->GetEntriesFast(); j++)
	{
		secTr = (EdbTrackP*)(secArr->At(j));
		if (priTr == secTr) {return 0; break;}
	}
	return priTr;
}
