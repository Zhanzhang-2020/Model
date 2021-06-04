
foldername = newArray(
	"ZZ_3D_tumor_0516_Jurkat_6h_dark"
	//"ZZ_3D_tumor_0501_Jurkat_0h"
	//"002", "003"
	);

raw_path = "./3D_tumoroid/"
date = '20210516'

for(ii=0; ii<lengthOf(foldername);ii++) {
	listdir = getFileList(raw_path + foldername[ii] );
	listdir = Array.sort(listdir);
	tifn = lengthOf(listdir);
	print("Porcessing "+foldername[ii]);
	for (iii=61;iii < 81;iii++) {
		G1 = raw_path+ foldername[ii]+ "/" + "well" + iii + "xy1" +"-c3.tiff";
		R1 = raw_path+ foldername[ii]+ "/" + "well" + iii + "xy1" +"-c4.tiff";
		//print(G1);
		open(G1);
		selectWindow("well" + iii + "xy1" +"-c3.tiff");
		run("Find focused slices", "select=55 variance=0.000 edge select_only verbose log");
		selectWindow("Find_Focus");
		a = getResult("Slice",0);
		b = getResult("Slice",nResults-1); //target z positions
		setAutoThreshold("Otsu dark");
		run("Threshold...");
		run("Create Selection");
		roi_name = iii;
		Roi.setName(roi_name);
		roiManager("Add");
		roiManager("Measure");
		selectWindow("Results");
		saveAs("Results", raw_path+"results"+"/"+date+"/"+foldername[ii] +"_area_"+iii+".csv");
		close("Results");
		//Roi.getBounds(xx, yy, width, height);
		open(R1);
		selectWindow("well" + iii + "xy1" +"-c4.tiff");
		run("Z Project...", "start="+ a +" stop=" + b +" projection=[Max Intensity]");
		roiManager("select", 0)
		//makeRectangle(xx, yy, width, height); //produce ROI based on previous green threshold
		run("Detect Particles", "ch1i ch1l ch1a=5 ch1s=5 rois=Ovals add=Nothing summary=Reset");
		close('Find_Focus');
		close('Threshold');
		close('Results');
		close('Log');
		//how to save the results window with only one row
		close("*");
		selectWindow("Summary");
		saveAs("Results", raw_path+"results"+"/"+date+"/"+ foldername[ii] +"_count_"+iii+".csv");
		close("Summary_"+iii+".csv");
		run("Select None");//Remove all the selections
		roiManager("reset");//Reset the ROI manager
		}
};
