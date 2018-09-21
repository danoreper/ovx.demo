Data for Schoenrock, 2016, "Ovariectomy results in inbred strain-specific increases in anxiety-like behavior in mice."

FST.csv:
Forced Swim Test (FST) Data. Each record represents FST data for a single mouse.

Fields:
-BATCH: the experimental batch.
-ANIMAL.ID: the unique integer ID of the mouse.
-STRAIN: The inbred strain of the mouse, one of 37 different strings.
-TREATMENT: Either OVX (ovariectomy) or SHAM (Sham surgery)
-ImbTime: In seconds, time that the mouse was immobile
-MobileTime: In seconds, time that the mouse was mobile.
-TotalTime: In seconds, time that the mouse was in the water.
-PctImb: Percent of the TotalTime that the mouse was immobile.


OF.csv:
Open field (OF) test data. Each record repreesnts OF data for a single mouse.

Fields:
-BATCH: the experimental batch.
-ANIMAL.ID: the unique integer ID of the mouse.
-STRAIN: The inbred strain of the mouse, one of 37 different strings.
-TREATMENT: Either OVX (ovariectomy) or SHAM (Sham surgery).
-TotDist: In centimeters, total distance travelled.
-MovTime: In seconds, time spent moving in the arena. 
-VMovNo: number of vertical movements
-MrgDist: in cm, distance travelled in the margin of the arena.
-CtrDist: In cm, distance travelled in the center of the arena.
-PctCenter: percent time spent in the center. 

