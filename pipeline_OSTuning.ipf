#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "Sarfia"

// Thomas Ryan 08/07/19 - Preprocessing step for raw data prior to ROI/segmentation by SARFIA or suite2p.
// All tif files in the given folder are loaded into Igor, channels are split and the imaging channel registered
// before saving the experiment and realigned tif with "_reg" suffix. Timewaves and locomotion waves are automatically
// generated and saved in the experiment file. Requires a few scripts from SARFIA
// Once this has completed, create the ROI maps for each registered movie, and import them to here (manually for now)
// then run the script postProcessOS()

// To Do - 	automate naming of experiment based on folder structure (line 57)
//				automate import of ROI masks from matlab data (in postprocess)
				// FINISH FITTING PROCEDURE
  
function preProcessOS()

	////////////////////////////////// some options \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	// Folder where your raw tif files are; use full path here.
	string pathstr="/Volumes/LocalDataHD/cmp34/Dropbox/Reg_Movies/Thiseas_022/24_07_19"

	// Whether you want to save the tif seperately (e.g. for suite2P analysis) or not
	variable saveTiff=1 	// use 0 if you do NOT want to save the tif. The image will still be stored in the experiment file (dependent on next option) and the raw data will be left alone

	// Whether you want to save the raw image files within the experiment. Significantly slims dow the experiment file and ou probably don't need it in the experiment if you are saving the tif seperately. Raw data is never touched
	variable saveExpImg=0 // save raw images
	variable saveCh1ExpImg=1 // save registered ch1
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	newpath/o path1, pathstr
	string filelist=indexedfile(path1,-1,".tif") 	// this will make a list of all tifs in the folder (careful it's only raw data)
	variable numFiles=itemsinList(filelist)		// simply counts how many files we need to loop through

	string message="Processing all tifs in folder " + pathstr
	print message

	// prepare some variable, strings and waves for the loop
	variable fsloop, nChannels, rL, RnChannels, dimZ, i, mm
	string fileN, regFileN, waveN, c1, c2, c3, regc1, stimName, locoName, expName

	for(fsloop=0;fsloop<numFiles;fsloop+=1)

		fileN=stringfromlist(fsloop,filelist)
		rL=strlen(fileN)

		waveN=fileN[0,(rL-5)]	// this is what the raw wave that we are loading will be called

		// internal wave names
		c1=waveN+"_ch1"			
		c2=waveN+"_ch2"
		c3=waveN+"_ch3"
		regc1=c1+"_reg"
		stimName="timewave_"+waveN
		locoName="loco_"+waveN
		
		// external file names 
		regFileN=regc1+".tif"		//We'll use this file name to make the tiff file later
		expName="expName.pxp"		
				
		// load the raw file
		message="Loading image file " + fileN + "; This is file " + num2str(fsloop+1) + " \ " + num2str(numFiles)
		print message
		AutoLoadScanImageSingle(pathstr, fileN)
		
		// reapply header information (this should make the deltas correct for all dimensions)
		duplicate/free/o $waveN Wave3D
		ApplyHeaderInfo(Wave3D)
		duplicate/o Wave3D $waveN
		
		// find the number of channels (only 3 channels compatible right now)
		duplicate/free/o/r=[][][0] $waveN firstrecpic
		redimension/n=(-1,-1) firstrecpic
		RnChannels =  nChannelsFromHeader(firstrecpic)

		// Split the channels and delete the raw movie (original .tifs are not touched)
		print "Splitting channels of recording"
		SplitChannels($waveN,RnChannels) 

		// Once split, there's no need to keep the (huge) unsplit imported data, so we'll clear it here.
		killwaves $waveN

		// Register Ch1; we can do this ourselves but currently using SARFIA script
		duplicate/free/o $c1 picwave
		RegisterStack(picwave, target=regc1)

		// Now make the stimulus timewave (Ch2); done by taking an average of each time frame. This will be binarized later in Leon's scripts so this is fine
		dimZ=dimsize($c2,2)							// count the time frames
		make/free/o/n=(dimZ) stimTrace			// make a new wave to store the timewave

		for(i=0;i<dimZ;i+=1)

			duplicate/o/free/r=[][][i] $c2 w
			wavestats/q w
			stimTrace[i]=V_avg

		endfor

		// subtract minimum to bring baseline to 0
		mm=wavemin(stimTrace)
		stimTrace-=mm

		// Now make the locomotion wave (Ch3) exactly the same way
		dimZ=dimsize($c3,2)				  		 	// count the time frames
		make/free/o/n=(dimZ) locoTrace			// make a new wave to store the timewave

		for(i=0;i<dimZ;i+=1)

			duplicate/o/free/r=[][][i] $c3 w
			wavestats/q w
			locoTrace[i]=V_avg

		endfor
		
		// subtract minimum to bring baseline to 0
		mm=wavemin(locoTrace)
		locoTrace-=mm

		// store the stim and loco waves in appopriately named waves
		duplicate/o stimTrace $stimName
		duplicate/o locoTrace $locoName

		// if saveTiff is set to 1, save the registered tiff in the same folder
		if(saveTiff!=0)
			duplicate/o/free $regc1 w
			imagesave/t="tiff"/s/ds=32/p=path1 w as regFileN
		
			// Now check if we want to keep the data in the experiment before we save it. 
			// N.B.Image data will NEVER be deleted if it is not saved as a tiff first (i.e. saveTiff=1).
			// BUT the other two channels WILL be deleted if the option saveExpImg is set to 0, 
			// EVEN IF saveTiff is set to 0
		
			if(saveCh1ExpImg!=1)
				killwaves $regc1
			endif
		
		endif
		
		if(saveExpImg!=1)
			killwaves $c1, $c2, $c3		
		endif
		
	endfor
	
	saveExperiment/p=path1 as expName
	
end

// Once ROI masks are imported, run this script to process all the movies in the experiment for OS analysis
// it will automatically apply the newest DFF0 and background subtraction procedures and perform the OSI analysis
// Once this is done, check that the analysed folders have suffixes of... 
//			"_1" if it was a 20 degree stimulus
//			"_2" if it was fullfield
 
function postProcessOS(nOri)

	variable nOri

	///////////////////OPTION\\\\\\\\\\\\\\\\\\\\\\
	variable deleteImages=1		// SET TO 0 IF YOU WANT TO KEEP ALL REGISTERED MOVIES IN THE EXPERIMENT!!
	////////////////////////////////////////////////	

	string regList=wavelist("*_reg",";","")		// find all the registered movies
	string ROIList=wavelist("RoiWave*",";","")		// find all the ROIWaves

	string reg, ROIn										// declare for individual files

	variable num=itemsinlist(regList), i			// count the number of movies

	for(i=0;i<num;i+=1)									// loop through each movie

		reg=stringfromlist(i,regList)					// find the name of the 'i'th registered movie...
		ROIn=stringfromlist(i,ROIList)					// ... and ROI wave

		DFF0($ROIn, $reg)									// Perform DFF0 and background subtraction

	endfor

	batchProcessOrientTuning(nOri)						// Process all movies for OSI analysis
	
	// clean up	
	wave destFrame,ROI,ExcluderegMovie,Histo,modAvCurve,tempFrame,tempFrameMovie,tempMovie,OriPref_All,OriPref_Pos_Hist
	killwaves destFrame,ExcluderegMovie,Histo,ROI,tempFrame,tempFrameMovie,tempMovie,OriPref_All,OriPref_Pos_Hist
	
	if(deleteImages!=0)
		string delList=wavelist("*_Ch1_reg",";","")
		killallwaveslist(delList)
	endif
end

// Thomas Ryan 14/04/19
// Computes the OSI of all ROIs within a field of view. Also correlates the ROI activity trace with the stimulus wave and stores
// the OSI values of positively correlated ROIs seperately. WARNING, USING THE BATCH SCRIPT WILL OVERWRITE THE ROI FOLDERS
// Inputs are simply the population wave following neuropil subtraction and DFF conversion.
// Outputs are organised by a seperate datafolder for each ROI and an "Analysis" folder for summary stats and tuning curves
// each containing:

// ... root:ROI_x:
// -> a matrix of the responses to repeats of each stimulus direction

// ... root:Analysis_[expName]
// -> results of CorrelationLinear procedure run against the stimulus including...
// * All correlation coefficients
// * All p-values
// * Identity of negatively ("[expName]_Corr_neg"), positively ("[expName]_Corr_pos"), and not ("[expName]_Corr_NS") correlated ROIs

// Pop-waves:
// -> the average of the integral across trials to each stimulus direction ("OrientTuningAve_Pop")
// -> the standard deviation of the integral across trials to each stimulus direction ("OrientTuningSD_Pop")
// -> OSI computed by Vector Average analysis of all ("OSI_VecAve_All") and only positively correlated ("OSI_VecAve_Pos") ROIs

// The stimulus, pop-wave and orientation order used to run the procedure is also stored in this folder
// OSI preferences and magnitudes can be computed by either fitting or Vector Average.
// DO NOT COMPARE/COMBINE HISTOGRAMS USING DIFFERENT METRICS 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	TO DO:New experiments will be recorded with 12 orientations and we may be able to fit; here are the orientations: 0,180,22.5,90,157.5,270,45,225,112.5,135,67.5,315 \\
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////BATCH FUNCTION\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function batchProcessOrientTuning(nOri)

	variable nOri
	variable method

	if(nOri>12)
		method=1
	else
		method =2
	endif

	string listDFF=waveList("*ch1_reg_DFF0",";","")
	string listStimW=wavelist("timewave*",";","")

	liststimw=sortlist(liststimw,";")
	listdff=sortlist(listdff,";")
	
	variable numList=itemsinlist(listDFF),i

	for(i=0;i<numList;i+=1)

		string NorPopName=stringfromlist(i,listDFF)
		string stimWaveName=stringfromlist(i,listStimW)

		OrientTuning(NorPopName,stimWaveName,method=method)

	endfor

end

////////////////MAIN FUNCTION\\\\\\\\\\\\\\\\\\\\\\\\\
function OrientTuning(NorPopName,stimWaveName[,method])

	variable method
	string NorPopName,stimWaveName
	
	duplicate/o/free $NorPopName NorPop
	
	if(paramisdefault(method))
	
		string mm = "What method would you like to use to compute OSI?"
		string METHODSTR="Vector Average (requires fewer samples in orientation space but not as accurate)"
		Prompt METHODSTR, mm, popup "Fitting (requires at least 12 orientations and 4 repeats);Vector Average (requires fewer samples in orientation space but not as accurate)"
		DoPrompt mm, METHODSTR
		
		if (stringmatch(METHODSTR,"Fitting (requires at least 12 orientations and 4 repeats)"))
			METHOD=1
		elseif (stringmatch(METHODSTR,"Vector Average (requires fewer samples in orientation space but not as accurate)"))
			METHOD=2
		endif
		
	endif
	
	// USER INPUT////////////////////////////////////////////////
	variable nReps=4
	string histograph = "yes" // Or "Yes"
	variable OSIThresh=0.4
	
	if(method==2)
		make/o/free orientations={0, 180, 90, 270, 45, 225, 135, 315}
	else
		make/o/free orientations={0,202.5,180,22.5,292.5,90,157.5,270,45,247.5,337.5,225,112.5,137.5,67.5,315}
	endif
	/////////////////////////////////////////////////////////////
	
	variable nFrames=dimsize(NorPop,0)
	variable nRois=dimsize(NorPop,1)
	variable deltat=dimdelta(NorPop,0)

	variable nOri=numpnts(orientations)
	make/o/free/n=(nOri) OriResp=0
	variable nStim=nOri*nReps
	
	duplicate/o/free $stimWaveName wStim


	setscale/p x,0,deltat,wStim
	
	// Before we do anything, check correlation with the stimulus wave to collect only those ROIs which are significantly correlated
	Correlation_Linear(stimwavename,NorPopName)
	
	// This creates a new datafolder called "Analysis" we will need a path to access later
	variable ll=strlen(NorPopName)
	string temp=	NorPopName[0,(ll-14)]
	
	string AnPath="root:" + temp
	
	variable m=mean(wStim),i,j,k,l

	// Find the start and end frames of each stimulation
	findlevels/Edge=1/Q/D=StimOn wStim, m // find where each stim comes on in seconds
	findlevels/Edge=2/Q/D=StimOff wStim, m // find where each stim goes off in seconds

	// Sort the orientations list so we can store the responses to each in the right order!
	duplicate/o/free orientations OriSorted
	sort OriSorted,OriSorted
	variable OriInt=OriSorted[1]-OriSorted[0]			// assumes distance between angles is consistent
	

	// Now extract the responses to each repeat of each orientation
	make/o/free/n=(nOri) OriResp=0	// free wave to store single response curves 
	make/o/free/n=(nOri,nReps) OriCurvesInt=0	// free wave to store all integral curves over repeats
	make/o/free/n=(nOri,nRois) OriAveCurve=0, OriStDevCurve	=0// These waves will store the average and SD of the integral waves across repeats
		
	// First, pull the trace and extract the segments for each stimulus repeat
	for(i=0;i<nRois;i+=1)
		duplicate/o/r=[][i] NorPop tempTrace	// pull full ROI trace; tempTrace
		redimension/n=(-1) tempTrace
		setscale/p x,0,deltat,tempTrace
		
		string RoiPath=AnPath+":ROI_" + num2str(i)	// create datafolder to store info aout this ROI
		NewDataFolder $RoiPath
		
		for(j=0;j<nOri;j+=1)	// cycle through the tested orientations

			variable thisOrient=orientations[j] // pull the actual orientation here
			string OriPath=RoiPath + ":Ori_" + num2str(thisOrient)	// create a wave by this name identifying the ROI and the Orientation
			
			for(k=0;k<nReps;k+=1)	// cycle through repeats

				variable thisStim=j+(k*nOri)	// identify which stimuli correspond to this orientation for repeat k
				variable onTime=StimOn[thisStim]	// identify the time this stimulus came on
				variable offTime=StimOff[thisStim]	// and off

				string outname= OriPath + "Reps"

				duplicate/o/r=(onTime,offTime) tempTrace stimTemp	// extract the time range for this stim from the trace
				setscale/p x,0,deltat,stimTemp
				
				if(k==0)
				
					make/o/free/n=(dimsize(stimTemp,0),nReps) repeatsTemp
					
				endif
			
				// check dimensions match; sometimes using time () instead of points [] can result in 10 or 11 frames extracted
					
				if((dimsize(repeatsTemp,0))!=(dimsize(stimTemp,0)))
					
					variable xx=x2pnt(NorPop,onTime)
					variable xxx=xx+(dimsize(repeatsTemp,0))
					duplicate/o/r=[xx,xxx] tempTrace stimTemp
					setscale/p x,0,deltat,stimTemp
				
				endif
					
				// Populate repeats matrix with all responses to all repeats to a given orientation in named waves
					
				repeatsTemp[][k]=stimTemp[p] 
				setscale/p x,0,deltat,repeatsTemp
				duplicate/o repeatsTemp $outName
				
				// Set negative values in the response to zero; i.e. only consider intensity values above 0 to count as the response. 
				// Hyperpolarisations and a drop to negative fluorescence may happen, but doesn't cancel out the response above 0
				
				variable npts=numpnts(stimtemp)
				
				for(l=0;l<npts;l+=1)
				
					if(stimtemp[l]<0)
					
						stimtemp[l]=0
						
					endif
					
				endfor
				
				// Compute the area under each curve for each repeat
				
				OriCurvesInt[j][k]=area(stimTemp)	// j is orientation (x) and k is repeats
				
				// now this OriCurvesInt wave contains the integral tuning curves for all repeats for each orientation
				
			endfor
		endfor
		
		// Still inside the ROI loop
		// Compute the mean and SD of the integral tuning curves into OriAveCurve and OriStDevCurve for this ROI
		
		matrixtranspose OriCurvesInt
		matrixOP/o/free Av=averagecols(OriCurvesInt)
		matrixtranspose OriCurvesInt
		matrixtranspose Av
		
		// Create a wave to store SD of the integral
		duplicate/o/free Av StDev
		StDev=0
			
		for(l=0;l<nOri;l+=1)
			duplicate/o/free/r=[l][] OriCurvesInt outTempCols
			wavestats/q outTempCols
			StDev[l]=V_Sdev
		endfor
			
		// Sort the waves so orientations are in correct order
		sort orientations,Av
		sort orientations,StDev
		
		//set orientations values correctly
		setscale/p x,0,OriInt, Av
		setscale/p x,0,OriInt, StDev
	
		OriAveCurve[][i]= Av[p]
		OriStDevCurve[][i]= StDev[p]
			
	endfor
	
	// store stuff in named waves
		
	string AveIntPopWave=AnPath + ":OrientTuningAve_Pop"
	string SDIntPopWave=AnPath + ":OrientTuningSD_Pop"	
	string OrientationsList=AnPath + ":Orientations"
	
	setscale/p x,0,OriInt, OriAveCurve
	setscale/p x,0,OriInt, OriStDevCurve
	
	duplicate/o OriAveCurve $AveIntPopWave
	duplicate/o OriStDevCurve $SDIntPopWave
	duplicate/o orientations $OrientationsList
		
	///////////////////////////////Raw////////////////////////////////
	// Create an average response (and their SD) to this stimulus over the repeats: 
	// these are secondary and currently stored in each ROIs individual folder in case we want them later
	string AveName=OriPath + "_Ave"
	string StDevName=OriPath + "_StDev"
						
	matrixTranspose stimTemp
	matrixOP/o/free Av=averagecols(stimTemp)
			
	// Create a wave to store SD
	duplicate/o Av StDev
	StDev=0
			
	for(l=0;l<numpnts(Av);l+=1)
		duplicate/o/free/r=[l][] repeatsTemp outTempCols
		wavestats/q outTempCols
		StDev[l]=V_Sdev
	endfor
			
	matrixtranspose StDev
	duplicate/o StDev $StDevName
	matrixTranspose Av
	duplicate/o Av $AveName
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Now we need to work out the OSI of each ROI
	// We can do this one of three ways; the Swindale method (Swindale, 1998) is the most accurate estimate but requires 
	// sampling orientation space with high resolution (at least 12 directions). fitting gaussians to the response curves
	// ("empirical OSI" a la Swindale, Shapley, Carandini and many many classical experiments) here generally results in a poor fit, 
	// most likely because we are not sampling enough of orientation space. Given that when it comes to fitting, GIGO (garbage in, 
	// garbage out), better to use a different measure. Fit-free methods include vector-sum or average and circular variance. 
	// Circular variance underestimates the OSI and is very sensitive to the width of the tuning curve. Vector analysis tends to overestimate,
	// but seems a good compromise with the number of repeats we have currently. 
	
	// To do: compute the circular variance as well...?
	// 
		 
	// Let's try to fit first... edit TR 27/04/19 - commented out preliminary fitting steps as clear that not enough orientations or repeats (or both) to fit properly
		 
	
	variable oriOrth
	make/o/n=(nRois) OSI_VecMag_All, OriPref_All
	OSI_VecMag_All=0		// list of OSI of all ROIs, regardless of response
	make/o/n=(200,nRois) vonMisesFit	// wave to store all fits

	duplicate/o/free/r=[][0] OriAveCurve tempAvCurve
	npts=(numpnts(tempAvCurve))/2
	make/o/n=(npts,nRois) modAvCurve_pop
	setscale/i x,0,180, modAvCurve_pop
	
	for(i=0;i<nRois;i+=1)
	
		// First, grab the curve and SD for this ROI
		duplicate/o/free/r=[][i] OriAveCurve tempAvCurve
		
		//duplicate/o/free/r=[][i] OriStDevCurve tempSDCurve
		redimension/n=(-1) tempAvCurve //tempSDCurve 

		// create a negative SD and Positive SD tuning curve for graphing
	
		//matrixop/o/free posSD=tempAvCurve+tempSDCurve
		//matrixop/o/free negSD=tempAvCurve-tempSDCurve

		//setscale/p x,0,OriInt, negSD
		//setscale/p x,0,OriInt, posSD
		
		// wrap the curve to 2*pi (so that responses to opposite directions are summed)
		
		npts=(numpnts(tempAvCurve))/2
		duplicate/o/free/r=[0,npts-1] tempAvCurve FirstmodAvCurve
		duplicate/o/free/r=[npts,npts*2] tempAvCurve SecondmodAvCurve
		matrixop/o modAvCurve=(FirstmodAvCurve+SecondmodAvCurve)/2
		
		// store the tuning curve
		
	
		variable OriIntRad=OriInt
		setscale/p x,0,OriInt, modAvCurve
		
		
		
		switch (METHOD)
		
			case 2:
				// Vector Average
				// grab the angles
				duplicate/o/r=[0,(nOri/2)-1] OriSorted theta
		
				// convert angles to radians
		
				theta[]= theta[p] * (pi/180)
				
				vectorAve(modAvCurve,theta)
				wave oriOut
				OriPref_All[i]=oriOut[0]
				OSI_VecMag_All[i]=oriOut[1]
				killwaves oriOut
				break
				
			case 1:
				// Fit vonMises curves to the tuning curves to find the preferred orientation
				setscale/i x,0,pi,modAvCurve	// convert x axis to radians
				Make/D/N=3/O W_coef	// estimate coefficient starting points. [0]=peak; [1]=preferred angle; [2]= width parameter
				wavestats/q modAvCurve			
				W_coef[0]=v_max
				w_coef[1]=pi/2
				w_coef[2]=v_maxloc
				
				FuncFit/q vonMises W_coef modAvCurve /D
				wave fit_modAvCurve
				//store curves and modded Av Curves
				vonMisesFit[][i]=fit_modAvCurve[p]
				modAvCurve_pop[][i]=modAvCurve[p]
				
				if(i==15)
					print w_coef[2]
				endif
				// find the peak and the peak x position (Ori_Pref) from the coefficients 
				OriPref_All[i]=w_coef[2]*(180/pi)
				
				if(OriPref_All[i]<0)
					OriPref_All[i]+=180
				endif
				
				//	find response at preferred ori
				variable prefResp=w_coef[0]
				
				// find orthogonal orientation by subtracting or adding 90 degrees
				variable anglePref=w_coef[1]/(pi/180)
				if(anglePref>=90)
					oriOrth=anglePref-90
				else
					oriOrth=anglePref+90
				endif
				variable div=pi/180
				oriOrth/=div
				// find response of model at oriOrth
				variable mmm=x2pnt(fit_modAvCurve,oriOrth)
				variable orthResp=fit_modAvCurve[mmm]
				
				//	OSI by measurement of the curve (PREFERRED_RESPONSE - ORTHOGONAL_RESPONSE)/PREFERRED_RESPONSE)
				OSI_VecMag_All[i]=(prefResp-orthResp)/prefResp
				
		endswitch
		
		if(method==1)
			string vonMisesName=AnPath + ":" + "vonMisesFit"
			setscale/i x,0,180,vonMisesFit
			duplicate/o vonMisesFit $vonMisesName
			killwaves root:vonMisesFit
			string modRespCurveName=AnPath + ":" + "modAvCurve_pop"
			duplicate/o modAvCurve_pop $modRespCurveName
			killwaves root:modAvCurve_pop
		endif
	endfor
	
	// The above vector average computation (method==2) does not consider a level of noise that could lead to high OSI values if the response is mostly noise. 
	//	Now extract only those ROIs which showed a significant positive correlation to the stimulus (i.e. show some response so we're not measuring noise)
	string NegCorrRoisName=AnPath + ":" + NorPopName + "_Corr_neg"
	duplicate/o/free $NegCorrRoisName NegROIlist
	variable numNegRois=numpnts(NegROIlist)
	
	string NSCorrRoisName=AnPath + ":" + NorPopName + "_Corr_NS"
	duplicate/o/free $NSCorrRoisName NSROIlist
	variable numNSRois=numpnts(NSROIlist)
	
	string PosCorrRoisName=AnPath + ":" + NorPopName + "_Corr_pos"
	duplicate/o/free $PosCorrRoisName ROIlist
	
	variable numPosRois=dimsize(ROIList,0)
	make/free/o/n=(numPosRois) OSI_VecMag_Pos, OriPref_Pos
	
	for(i=0;i<numPosRois;i+=1)
	
		variable tempMag=OSI_VecMag_All[(RoiList[i])]
		OSI_VecMag_Pos[i]=tempMag
		
		variable tempPref=OriPref_All[(RoiList[i])]
		OriPref_Pos[i]=tempPref
		
	endfor
	
	// Move results waves to analysis folder - Magnitude

	string newname=AnPath + ":OSI_VecMag_All"
	duplicate/o OSI_VecMag_All $newName
	
	newname=AnPath + ":OSI_VecMag_Pos"
	duplicate/o OSI_VecMag_Pos $newName
	
	if(stringmatch(histograph,"Yes"))
		// Make the histogram of vector magnitudes
		
		Make/N=25/O OSI_VecMag_Pos_Hist;DelayUpdate
		Histogram/B={0,0.05,20} $newName,OSI_VecMag_Pos_Hist;DelayUpdate
		newName=AnPath + ":OSI_VecMag_Pos_Hist"
		duplicate/o OSI_VecMag_Pos_Hist $newName
		Display $newName
		ModifyGraph mode=5,hbFill=4,rgb=(0,0,0)
		Label left "\\Z16No. ROIs";DelayUpdate
		Label bottom "\\Z16OSI"
	endif
	
	
	// Now of the Orientation preferences

	newname=AnPath + ":OriPref_All"
	duplicate/o OriPref_All $newName
	
	newname=AnPath + ":OriPref_Pos"
	duplicate/o OriPref_Pos $newName
	if(stringmatch(histograph,"Yes"))
		Make/N=36/O OriPref_Pos_Hist;DelayUpdate
		Histogram/B={0,10,18} $newName,OriPref_Pos_Hist;DelayUpdate
		newName=AnPath + ":OriPref_Pos_Hist"
		duplicate/o OriPref_Pos_Hist $newName
		Display $newName
		ModifyGraph mode=5,hbFill=4,rgb=(0,0,0)
		Label left "\\Z16No. ROIs";DelayUpdate
		Label bottom "\\Z16Orientation angle (degrees)"
	
	endif
	
	// Print useful metrics
	
	// Number of negative cells 
	string message= "There are " + num2str(numNegRois) + " / " + num2str(nRois) + " (" + num2str((numNegRois/nRois)*100) + "%) ROIs negatively correlated with the stimulus"
	print message

	// not correlated/no activity
	message= "There are " + num2str(numNSRois) + " / " + num2str(nRois) + " (" + num2str((numNSRois/nRois)*100) + "%) ROIs showed no correlation or no activity"
	print message

	// Number of responsive cells (correlated with stimulus)
	message= "There are " + num2str(numPosRois) + " / " + num2str(nRois) + " (" + num2str((numPosRois/nRois)*100) + "%) ROIs positively correlated with the stimulus"
	print message
	print "Of these..."

	// How many have an OSI above threshold?
	variable numOSRois = 0
	for(i=0;i<numPosRois;i+=1)
		if(OSI_VecMag_Pos[i]>=OSIThresh)
			numOSRois+=1
		endif
	endfor
	message=num2str(numOSRois) + " of the ROIs responsive to the stimulus have an estimated OSI above a user defined threshold of " + num2str(OSIThresh)
	print message
	
	Make/O/N=(4,2) Numbers
	
	Numbers[0][0]=numNSRois
	Numbers[1][0]=numPosRois	
	Numbers[2][0]=numNegRois
	Numbers[3][0]=nRois		
	
	Numbers[0][1]=(numNSRois/nRois)*100
	Numbers[1][1]=(numPosRois/nRois)*100	
	Numbers[2][1]=(numNegRois/nRois)*100
	Numbers[3][1]=(nRois/nRois)*100		
	
	
	SetDimLabel 0, 0, NotSign, Numbers
	SetDimLabel 0, 1, Pos, Numbers
	SetDimLabel 0, 2, Neg, Numbers
	SetDimLabel 0, 3, All, Numbers    

	SetDimLabel 1, 0, Raw, Numbers        //Set column labels
	SetDimLabel 1, 1, Percent, Numbers
    
	newname=AnPath + ":Numbers"
	duplicate/o Numbers $newName
	//cleanup
	

	
	
	wave W_StatsLinearCorrelationTest
	killwaves OSI_VecMag_Pos_Hist,Stimoff,StimOn,stimTemp,tempTrace,theta,W_StatsLinearCorrelationTest, OSI_VecMag_Pos, w_coef, modAvCurve,fit_modAvCurve 
end

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Written by Rozan Vroman, modified Leon Lagnado 16/04/19, modified Thomas Ryan 29/04/19//
// Correlates each ROI from a pop wave with a stimulus wave to identify responsive neurons

Function Correlation_linear(LWN, MWN)
	string LWN, MWN

	duplicate/FREE/o $MWN All
	duplicate/FREE/o $LWN Stim
	string Analysis
	variable ll=strlen(MWN)
	string 	MWN_c=MWN[0,(ll-14)]
	Analysis="root:"+MWN_c
	newdatafolder $(Analysis)
	variable sigpvalue = 0.02
	variable nshuffles = 1000
	variable pvalue = 0
	variable i

	////// GUI ///////
	string newWaveName
	string SPEARorPEAR_select,histograph
	variable Hz=5

	histograph = "Yes please :)" // Or "No, thanks"
	Hz = 1/dimsize($MWN, 0)

	make/FREE/o/n=(dimsize(All,1)) List_Corr_STIM_CCvalues
	make/FREE/o/n=(dimsize(All,1)) List_Corr_STIM_Pvalues

	Make/FREE/o/n=(nshuffles) RandomWave 
	Make/FREE/o/n=(nshuffles) RandomCorValue

	Make/FREE/o/n=0 ROI_list_Corr_NS
	Make/FREE/o/n=0 ROI_list_Corr_pos
	Make/FREE/o/n=0 ROI_list_Corr_neg

	variable LenRandom=dimsize(RandomWave,0)
	variable LenAll=dimsize(All,0)
	variable deltat=dimdelta(All,0)
	variable LenROIs=dimsize(All,1)

	// Big Loop through all ROIs
	MatrixOp/O All_NANS=ReplaceNaNs(aLL,0)

	FOR (i=0;i<LenROIs;i+=1)
		duplicate/FREE/o/r=[][i] All_NANS temp1
		StatsLinearCorrelationTest/ALPH=0.05/Q stim, temp1
		wave W_StatsLinearCorrelationTest
		List_Corr_Stim_CCvalues[i]=W_StatsLinearCorrelationTest[1] // 1 for Pearson
		List_Corr_Stim_Pvalues[i]=W_StatsLinearCorrelationTest[7]  //p value	
	ENDFOR

	// Devide the ROIs in significantly positively or negatively correlated with locomotion, or not significant
	FOR (i=0;i<LenROIs;i+=1)
		if ((List_Corr_Stim_Pvalues[i] > sigpvalue) && (List_Corr_Stim_Pvalues[i] < (1-sigpvalue)))
			InsertPoints/M=0 (dimsize(ROI_list_Corr_NS,0)),1, ROI_list_Corr_NS
			ROI_list_Corr_NS[dimsize(ROI_list_Corr_NS,0)-1]=i
		elseif ((List_Corr_Stim_Pvalues[i] < sigpvalue) && List_Corr_Stim_CCvalues[i]>=0)
			InsertPoints/M=0 (dimsize(ROI_list_Corr_pos,0)),1, ROI_list_Corr_pos
			ROI_list_Corr_pos[dimsize(ROI_list_Corr_pos,0)-1]=i
		elseif ( (List_Corr_Stim_Pvalues[i]< sigpvalue) && List_Corr_Stim_CCvalues[i]<0)
			InsertPoints/M=0 (dimsize(ROI_list_Corr_neg,0)),1, ROI_list_Corr_neg
			ROI_list_Corr_neg[dimsize(ROI_list_Corr_neg,0)-1]=i
		endif
	ENDFOR


	string window_name2=winname(0,2,1)
	//killwindow $window_name2
	//killwaves W_StatsLinearCorrelationTest, tempRandom, Loco_Corrected

	// Put the new waves in the "Analysis" Folder
	newWaveName = Analysis + ":"+ MWN
	duplicate/o All $newWaveName

	newWaveName = Analysis + ":"+ LWN
	duplicate/o Stim $newWaveName

	newWaveName = Analysis + ":"+ MWN + "_Corr_neg"
	duplicate/o ROI_list_Corr_neg $newWaveName

	newWaveName = Analysis + ":"+ MWN + "_Corr_pos"
	duplicate/o ROI_list_Corr_pos $newWaveName

	newWaveName = Analysis + ":"+ MWN + "_Corr_NS"
	duplicate/o ROI_list_Corr_NS $newWaveName

	newWaveName = Analysis + ":"+ MWN + "_Corr_CC"
	duplicate/o List_Corr_Stim_CCvalues $newWaveName

	newWaveName = Analysis + ":"+ MWN + "_Corr_P"
	duplicate/o List_Corr_Stim_Pvalues $newWaveName
	
	//if (stringmatch(histograph,"Yes please :)"))
	//	CorrelationHistPlots(MWN)
	//endif
END

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Thomas Ryan 17/04/19 - Implementation of vector average quantification of orientation selectivity index for a given ROI
// Equation from Kerlin AM, Andermann ML, Berezovskii VK, Reid RC. Broadly tuned response properties of diverse inhibitory neuron subtypes in mouse visual cortex. Neuron. 2010;67(5):858–871. doi:10.1016/j.neuron.2010.08.002 
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3327881/#R78. Originally proposed by Swindale 1997 https://www.ncbi.nlm.nih.gov/pubmed/9518026

// Equation in note form: OSI = ((ΣR(θi)sin(2θi))^2 + (ΣR(θi)cos(2θi))^2)^1/2/ ΣR(θi), where θi is the orientation of each stimulus and R(θi) is the response to that stimulus.


function vectorAve(rVec,theta)

	wave rVec, theta
	
	//variable mVec=wavemin(rVec)
	//rVec-=mVec								// Not sure we should subtract this; it amplifies the OSI but assumes that there is a basal response to any stimulus that we do not consider in computing the selectivity. Not sure this is fair given that DF/F SHOULD reflect a difference from a true baseline outside of the stimulus periods. I guess it just depends on our definition of OSI
	variable nTheta= dimsize(rVec,0)
	
	variable i, prefOri
	
	duplicate/o/free theta, a, b
	
	for (i=0;i<nTheta;i+=1)
	
		variable dubTheta=2*theta[i]
		variable sindubTheta= sin(dubTheta)
		variable cosdubtheta=	cos(dubTheta)	
	
		b[i] = rVec[i] * sindubTheta
		a[i] = rVec[i] * cosdubTheta
		
	endfor
			
	variable sumb= sum(b)
	variable suma = sum(a)
	variable sumVec=sum(rVec)
	variable sumOfSquare=sqrt(suma^2 + sumb^2)
	variable vecMag = sumOfSquare/sumVec
	
	// preferred orientation can be found in radians as 0.5arctan(sumSin/sumCos) if sumCos>0
	// OR 180+0.5arctan(sumSin/sumCos) if sumCos<0
	variable div=sumb/suma
	
	prefOri = 0.5*atan(div)

	prefOri*=(180/pi)
	
	if(suma<0)
		prefOri+=90
	elseif(sumb<0)
		prefOri+=180
	endif
	
	make/o/n=2 oriOut
	oriOut[0]=prefOri
	oriOut[1]=vecMag
	
end

//////////////////////////////////////////////////////////////////////////////////////
// Thomas Ryan 15/04/19 - short function to display the full tuning curve of ROI n

function displayTuningCurve(expName,ROIn)
	string expName
	variable Roin
	variable rep=1 	// set to 	0- display graphs
	//				1- update graphs

	string aveCurve="root:"+expName+":modAvCurve_Pop"
	string OSIVecMag_All="root:"+expName+":OSI_VecMag_All"
	string OriPref_All="root:"+expName+":OriPref_All"
	duplicate/o/r=[][Roin] $aveCurve tempAvCurve

	// pull the fit
	
	string vonMises="root:"+expName+":vonMisesFit"
	duplicate/o/r=[][Roin] $vonMises tempFit
	
	if(rep==0)

		display tempAvCurve
		appendtograph tempFit
	endif
	duplicate/o/free $OSIVecMag_All OSIVecMag_Allwave
	duplicate/o/free $OriPref_All OriPref_Allwave
	
	variable tempOSI=OSIVecMag_Allwave[Roin]
	variable tempPref=OriPref_Allwave[Roin]
	print "OSI = " + num2str(tempOSI)
	print "Preferred Orientation = " + num2str(tempPref)
end


/////////////////////////////////////////////////////////
// Thomas Ryan 15/04/19 - short function to access the response curves of given ROI

function displayResponses(ROIn,Angle)

	variable ROIn,Angle

	variable rep=1 // set to 1 if cycling through ROIs to avoid, 0 to display new figures each run

	string path="root:ROI_" + (num2str(ROIn)) + ":"

	string repsN= path + "Ori_" + (num2str(Angle)) + "Reps"
	string AveN= path + "Ori_" + (num2str(Angle)) + "_Ave"
	string stDevN= path + "Ori_" + (num2str(Angle)) + "_StDev"

	duplicate/o $repsN tempRepsDisp
	duplicate/o $AveN tempAveDisp
	duplicate/o $stDevN tempStDevDisp

	matrixOP/o tempSDPosDisp=tempAveDisp+tempStDevDisp
	matrixOP/o tempSDNegDisp=tempAveDisp-tempStDevDisp

	variable numReps=dimsize(tempRepsDisp,1),i
	string thisRep

	for(i=0;i<numReps;i+=1)

		thisRep="tempRep_" + (num2str(i+1)) + "_Disp"
		duplicate/o/r=[][i] tempRepsDisp $thisRep

		if(i==0)
			display $thisRep
		else
			appendtograph $thisRep
		endif

	endfor

	appendtograph tempAveDisp

end



/////// Thomas Ryan 24/05/19 - this function pools all OSI data within a given 
// experiment that contains the info to make histograms for a given timepoint. Make sure that
// analysed folders have the "_1" suffix if the stimulus was 20degree, and _2 if it was fullfield. 

function poolDataWithin()

	variable hist=2 // 2 to produce histograms and display them, 1 to produce histograms but not display them, 0 for neither

	variable numConditions=2,i,j
	string dat=dataFolderDir(1) // get names of all folders in the experiment; this produced a comma seperated list with a semicolon at he end (which we need to remove)

	dat=removeending(dat,";")

	string listCondition_1=ListMatch(dat,"*_1",",")	// pull folders which end in _1, specifying a common condition
	string listCondition_2=ListMatch(dat,"*_2",",")

	variable numFolders_1=itemsinlist(listCondition_1,",")
	variable numFolders_2=itemsinlist(listCondition_2,",")
	
	string list_OriPrefAll="", list_OriPrefPos="", list_OSIMagPos="", list_OSIMagAll=""
	
	for(i=0;i<numFolders_1;i+=1)
	
		string folderName=stringfromlist(i,listCondition_1,",") // find the folder name
		string rFolderName="root:"+folderName +":"	// add it to the root folder name
		string path_OriPrefAll= rFolderName + "OriPref_All" // ID lists of data to concatenate
		string path_OriPrefPos= rFolderName + "OriPref_Pos"
		string path_OSIMagPos= rFolderName + "OSI_VecMag_Pos" 
		string path_OSIMagAll= rFolderName + "OSI_VecMag_All"
	
		string target_OriPrefAll= "OriPrefAll_" + (num2str(i))	// define targets in the root folder for these waves
		string target_OriPrefPos= "OriPref_Pos_" + (num2str(i))
		string target_OSIMagPos= "OSI_VecMag_Pos_" + (num2str(i))
		string target_OSIMagAll= "OSI_VecMag_All_" + (num2str(i))
	
		duplicate/o $path_OriPrefAll $target_OriPrefAll	// duplicate source wave from each folder to target
		duplicate/o $path_OriPrefPos $target_OriPrefPos
		duplicate/o $path_OSIMagPos $target_OSIMagPos
		duplicate/o $path_OSIMagAll $target_OSIMagAll

		list_OriPrefAll=addlistitem(	target_OriPrefAll, list_OriPrefAll,";")	// create a list of each type of item to concatenate outside the loop
		list_OriPrefPos=addlistitem(	target_OriPrefPos, list_OriPrefPos,";")
		list_OSIMagPos=addlistitem(	target_OSIMagPos, list_OSIMagPos,";")
		list_OSIMagAll=addlistitem(	target_OSIMagAll, list_OSIMagAll,";")
		
	endfor
	
	print "All orientation preferences : " + list_OriPrefAll
	print ""
	print "Stim responsive orientation preferences : " + list_OriPrefPos
	print ""
	print "All OSI magnitudes : " + list_OSIMagAll
	print ""
	print "Stim responsive OSI magnitudes : " + list_OSIMagPos	
	print ""
	
	concatenate/kill/np/o list_OriPrefAll, concat_OriPrefAll_20deg 
	concatenate/kill/np/o list_OriPrefPos, concat_OriPrefPos_20deg
	concatenate/kill/np/o list_OSIMagPos, concat_OSIVecMagPos_20deg
	concatenate/kill/np/o list_OSIMagAll, concat_OSIVecMagAll_20deg
	
	list_OSIMagAll="",list_OriPrefPos="",list_OSIMagPos="",list_OriPrefAll=""	// reset lists
	
	
	for(i=0;i<numFolders_2;i+=1)
	
		folderName=stringfromlist(i,listCondition_2,",") // find the folder name
		rFolderName="root:"+folderName +":"	// add it to the root folder name
		path_OriPrefAll= rFolderName + "OriPref_All" // ID lists of data to concatenate
		path_OriPrefPos= rFolderName + "OriPref_Pos"
		path_OSIMagPos= rFolderName + "OSI_VecMag_Pos" 
		path_OSIMagAll= rFolderName + "OSI_VecMag_All"
	
		target_OriPrefAll= "OriPrefAll_" + (num2str(i))	// define targets in the root folder for these waves
		target_OriPrefPos= "OriPref_Pos_" + (num2str(i))
		target_OSIMagPos= "OSI_VecMag_Pos_" + (num2str(i))
		target_OSIMagAll= "OSI_VecMag_All_" + (num2str(i))
	
		duplicate/o $path_OriPrefAll $target_OriPrefAll	// duplicate source wave from each folder to target
		duplicate/o $path_OriPrefPos $target_OriPrefPos
		duplicate/o $path_OSIMagPos $target_OSIMagPos
		duplicate/o $path_OSIMagAll $target_OSIMagAll

		list_OriPrefAll=addlistitem(	target_OriPrefAll, list_OriPrefAll,";")	// create a list of each type of item to concatenate outside the loop
		list_OriPrefPos=addlistitem(	target_OriPrefPos, list_OriPrefPos,";")
		list_OSIMagPos=addlistitem(	target_OSIMagPos, list_OSIMagPos,";")
		list_OSIMagAll=addlistitem(	target_OSIMagAll, list_OSIMagAll,";")
	endfor
	
	print "All orientation preferences (FF) : " + list_OriPrefAll
	print ""
	print "Stim responsive orientation preferences (FF) : " + list_OriPrefPos
	print ""
	print "All OSI magnitudes (FF) : " + list_OSIMagAll
	print ""
	print "Stim responsive OSI magnitudes (FF) : " + list_OSIMagPos	
	print ""
	
	// concatenate all data from the lists

	concatenate/kill/np/o list_OriPrefAll, concat_OriPrefAll_FF 
	concatenate/kill/np/o list_OriPrefPos, concat_OriPrefPos_FF
	concatenate/kill/np/o list_OSIMagPos, concat_OSIVecMagPos_FF
	concatenate/kill/np/o list_OSIMagAll, concat_OSIVecMagAll_FF
		
	// check for negatives and add 180 to them if needed
	removeNegs("concat_OriPrefAll_20deg")
	removeNegs("concat_OriPrefPos_20deg")
	removeNegs("concat_OriPrefAll_FF")
	removeNegs("concat_OriPrefPos_FF")
	
	string	listOriHists="concat_OriPrefAll_20deg_hist;concat_OriPrefPos_20deg_hist;concat_OriPrefAll_FF_hist;concat_OriPrefPos_FF_hist;"
	string listMagHists="concat_OSIVecMagPos_20deg_hist;concat_OSIVecMagAll_20deg_hist;concat_OSIVecMagPos_FF_hist;concat_OSIVecMagAll_FF_hist;"
	
	variable numOriList=itemsinlist(listOrihists)
	variable numMagList=itemsinlist(listMagHists)
	/////// We now have 8 lists of numbers containing the orientation preferences and magnitudes of all experiments within each condition//
	// Next we want to make histograms of these (only ROIs positively correlated with stimulus for now)
	
	if(hist!=0)
		//make histograms
	
		Make/N=18/O concat_OriPrefAll_20deg_Hist,concat_OriPrefPos_20deg_Hist,concat_OriPrefAll_FF_Hist,concat_OriPrefPos_FF_Hist
		Make/N=20/O concat_OSIVecMagPos_20deg_Hist,concat_OSIVecMagAll_20deg_Hist,concat_OSIVecMagPos_FF_Hist,concat_OSIVecMagAll_FF_Hist
		
		Histogram/B={0,10,18} concat_OriPrefAll_20deg,concat_OriPrefAll_20deg_Hist;DelayUpdate
		Histogram/B={0,10,18} concat_OriPrefPos_20deg,concat_OriPrefPos_20deg_Hist;DelayUpdate
		Histogram/B={0,10,18} concat_OriPrefAll_FF,concat_OriPrefAll_FF_Hist;DelayUpdate
		Histogram/B={0,10,18} concat_OriPrefPos_FF,concat_OriPrefPos_FF_Hist;DelayUpdate
		
		Histogram/B={0,0.05,20} concat_OSIVecMagPos_20deg,concat_OSIVecMagPos_20deg_Hist;DelayUpdate
		Histogram/B={0,0.05,20} concat_OSIVecMagAll_20deg,concat_OSIVecMagAll_20deg_Hist;DelayUpdate
		Histogram/B={0,0.05,20} concat_OSIVecMagPos_FF,concat_OSIVecMagPos_FF_Hist;DelayUpdate
		Histogram/B={0,0.05,20} concat_OSIVecMagAll_FF,concat_OSIVecMagAll_FF_Hist;DelayUpdate

		if(hist==2)
			//display histograms
			for(j=0;j<numOriList;j+=1)
				string histstr=stringfromlist(j,listorihists)
				Display $histstr
				ModifyGraph mode=5,hbFill=5,rgb=(0,0,0)
				Label left "Frequency";DelayUpdate
				Label bottom "Orientation (mod 180) (degrees)"
			endfor
			
			for(j=0;j<numMagList;j+=1)
				histstr=stringfromlist(j,listmaghists)
				Display $histstr
				ModifyGraph mode=5,hbFill=5,rgb=(0,0,0)
				Label left "Frequency";DelayUpdate
				Label bottom "Orientation Selectivity Index"
			endfor
		endif
	endif
	
end


function removeNegs(vecName)
	string vecName
	
	duplicate/o/free $vecName wW
	variable num=numpnts(wW),i
	for(i=0;i<num;i+=1)
		if(wW[i]<0)
			wW[i]+=180
		endif
	endfor

	duplicate/o wW $vecName
end

function killAllWavesList(List)
	string List
	Variable numList = ItemsInList(list),j
			
	for (j=0;j<numList;j+=1)
		KillWaves/Z $(StringFromList(j, list))
	endfor
	
end


function vectorAve2(rVec,theta) // this is a version of this function that cuts the additional orientations for back compatability with previous experiments with only 8 orientations tested.

	wave rVec, theta
	
	//variable mVec=wavemin(rVec)
	//rVec-=mVec											// Not sure we should subtract this; it amplifies the OSI but assumes that there is a basal response to any stimulus that we do not consider in computing the selectivity. Not sure this is fair given that DF/F SHOULD reflect a difference from a true baseline outside of the stimulus periods. I guess it just depends on our definition of OSI
	
	make/o/n=4 rVecLeg, theta
	
	variable nTheta= dimsize(rVec,0)
	
	variable i, prefOri
	
	duplicate/o/free theta, a, b
	
	for (i=0;i<nTheta;i+=1)
	
		variable dubTheta=2*theta[i]
		variable sindubTheta= sin(dubTheta)
		variable cosdubtheta=	cos(dubTheta)	
	
		b[i] = rVec[i] * sindubTheta
		a[i] = rVec[i] * cosdubTheta
		
	endfor
			
	variable sumb= sum(b)
	variable suma = sum(a)
	variable sumVec=sum(rVec)
	variable sumOfSquare=sqrt(suma^2 + sumb^2)
	variable vecMag = sumOfSquare/sumVec
	
	// preferred orientation can be found in radians as 0.5arctan(sumSin/sumCos) if sumCos>0
	// OR 180+0.5arctan(sumSin/sumCos) if sumCos<0
	variable div=sumb/suma
	
	prefOri = 0.5*atan(div)

	prefOri*=(180/pi)
	
	if(suma<0)
		prefOri+=90
	elseif(sumb<0)
		prefOri+=180
	endif
	
	make/o/n=2 oriOut
	oriOut[0]=prefOri
	oriOut[1]=vecMag
	
end

//LL 15/4/2019
//Input is the ROI mask already generated from Suite2P and the registered movie.
//Does interpolation to make movie square.  Corrects background and calculates DF/F0.  
//Leaves behind pop wave with subscript "_DFF0" and the squared up movie ("_sqr") and corresponding ROI mask ("sqrROI").
//Approach is to measure local bkg and implement Svoboda factor.  This prevents -ve signals and is sensible.
//The local bkg is calculated form mask that approximates all ROIs as circles.  
//Measures average radius of ROIs and excludes pixels that are within 1.5*radius (i.e those that are too close) 
//but includes those that are > 1.5 r and < 2.25 r, as log as they are not on a neighbouring ROI exclusion zone.

//To do: tidy up edges of squared ROI mask that are not integers

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

Function DFF0(ROIMask, registeredMovie)

	wave ROIMask, registeredMovie

	wavestats/q ROIMask
	variable nROIs = (V_min*(-1))+1						//nROIs is the number of ROIs - edit - TR - 15/04/19 ROI index begins at 0: +1 to get nROI
	variable nx = dimsize(ROIMask, 0)
	variable ny = dimsize(ROIMask,1)
	variable nFrames = dimsize(registeredMovie,2)
	variable deltat = dimdelta(registeredMovie, 2)
	variable exclusionradius						//Distance from center of each ROI to be excluded. Currently in pixels
	variable maxradius								//Max distance from center of each ROI to be included. Currently in pixels
	variable SvobodaFactor = 0.5					//Otherwise, bkg seems to be overestimated. 
	variable scalefactor

	Make/O/N=(nFrames, nROIs) ROI_pop=0, BKG_pop=0
	Make/O/N=(nFrames) ROI
	SetScale /P x,0,deltat, ROI_pop, BKG_pop, ROI
	Make/O/N=(nROIs) ROInPixels, SurroundnPixels
	Make/O/N=(nROIs, 2)  ROIcentres=0

	String inputwavename = NameofWave(registeredMovie)
	String DFFOwavename = inputwavename + "_DFF0"					//This will be the output wave
	String Squarewavename = inputwavename + "_sqr"	
	String ROIMaskwavename = inputwavename + "_sqrROI"

	variable i, j, k, l, n, ROInum
	variable m = 0									//Pixel value of ROI number x will be -x
	variable F0 									//This will be the fluorescence at rest
	variable, dx, dy


	//Rescale movie to be square
	//Make a temp movie which is downscaled appropriately
	Make/O/N=(nx, ny) tempFrameMovie
	if(nx>ny)
		Make/O/N=(ny, ny, nFrames) tempMovie
		Make/O/N=(ny, ny) tempMask, tempFrame, destFrame
		scalefactor=ny/nx
	else
		Make/O/N=(nx, nx, nFrames) tempMovie
		Make/O/N=(nx, nx) tempMask, tempFrame, destFrame
		scalefactor=nx/ny
	endif

	for(i=0; i<nFrames; i+=1)
		tempFrameMovie[][]=registeredMovie[p][q][i]
		if(nx>ny)
			ImageInterpolate/F={scalefactor, 1}/DEST=destFrame Bilinear tempFrameMovie
		else
			ImageInterpolate/F={1, scalefactor}/DEST=destFrame Bilinear tempFrameMovie
		endif
		tempMovie[][][i]= destFrame[p][q]
	endfor

	Duplicate/O tempMovie regMovie

	//Rescale ROI mask to be square
	if(nx>ny)
		ImageInterpolate/F={scalefactor, 1}/DEST=destFrame Bilinear ROIMask
	else
		ImageInterpolate/F={1, scalefactor}/DEST=destFrame Bilinear ROIMask
	endif
	tempMask[][]= destFrame[p][q]
	
	Duplicate/O tempMask sqROIwave, surround_ROI
	surround_ROI = 1
	KillWaves/Z/F tempMask

	// insert TR - 15/04/19 - reset nx and ny variables to be the dimensions of the squared waves

	nx = dimsize(sqROIWave, 0)
	ny = dimsize(sqROIWave,1)

	//Find the center of mass of each ROI
	ROInPixels=0
	for (i=0; i< nx; i+=1)
		for (j=0; j< ny; j+=1)
			ROInum = sqROIwave[i][j] *-1
			if(ROInum >= 0)
				ROInPixels[ROInum]+=1	
				ROIcentres[ROInum][0]+=i
				ROIcentres[ROInum][1]+=j
			endif
		endfor
	endfor

	for (i=0;i<nROIs;i+=1)	
		ROIcentres[i][]/=ROInPixels[i]	
	endfor

	//Calculate the concentric circles over each ROI for the exclusion mask
	variable aveROIradius = sqrt(mean(ROInPixels)/Pi)			//mean radius in pixels
	exclusionradius = aveROIradius*1.25
	maxradius	 = exclusionradius*1.5

	//Make the exclusion mask
	Duplicate/O sqROIwave ExclusionMask
	ExclusionMask = 1

	for (i=0; i< nx; i+=1)
		for (j=0; j< ny; j+=1)
			for (k=0;k<nROIs;k+=1)
				dx = abs(i-ROIcentres[k][0])
				dy = abs(j-ROIcentres[k][1])
				if( (sqrt(dx^2 + dy^2) < exclusionradius))
					ExclusionMask[i][j]=0
				endif
			endfor
		endfor
	endfor

	//Make a version of the movie in which exclusion zones are blanked to zero
	// edit TR - 15/04/19 - Previously only multiplied the first frame... now does all frames in a loop

	Duplicate/O regMovie ExcluderegMovie

	for(i=0;i<(dimsize(ExcluderegMovie,2));i+=1)
		duplicate/o/free/r=[][][i] ExcluderegMovie tempExc
		tempExc*=ExclusionMask
		ExcluderegMovie[][][i]=tempExc[p][q]
	endfor

	//Now calculate signal in each ROI
	ROInPixels[]=0
	for (i=0; i< nx; i+=1)
		for (j=0; j< ny; j+=1)
			ROInum = sqROIwave[i][j]
			if(ROInum <= 0)
				ROInPixels[ROInum*-1]+=1
				MatrixOP/O ROI = beam(regMovie,i,j)
				ROI_pop[][(ROInum*-1)] += ROI[p]    
			endif
		endfor
	endfor

	for(i=0; i<nROIs; i+=1)
		ROI_pop[][i]= ROI_pop[p][i]/ROInPixels[i]
	endfor

	//Measure the surround
	SurroundnPixels[]=0
	for (i=0; i< nx; i+=1)
		for (j=0; j< ny; j+=1)
			if(ExclusionMask[i][j]==1)
				for (k=0;k<nROIs;k+=1)	
					dx = abs(i-ROIcentres[k][0])
					dy = abs(j-ROIcentres[k][1])
					if((sqrt(dx^2 + dy^2) < maxradius))
						SurroundnPixels[k]+=1
						MatrixOP/O ROI = beam(regMovie,i,j)   
						BKG_pop[][k] += ROI[p] 						 
					endif
				endfor
			endif
		endfor
	endfor

	BKG_pop*=SvobodaFactor
	for(i=0; i<nROIs; i+=1)
		BKG_pop[][i]= BKG_pop[p][i]/SurroundnPixels[i]
	endfor	

	//Now background correct (i.e subtract background)
	Duplicate/O ROI_pop BKGCorr_pop DFF0_pop
	matrixop/o BKGCorr_pop = ROI_pop - BKG_pop

	//Now find F0 and calculate DF/F0
	// edit TR - 15/04/19 - now computes f0 as the mode rather than the average
	for(i=0; i<nROIs; i+=1)
		duplicate/free/o/r=[][i] BKGCorr_pop ROI
		
		Make/o/n=100 Histo
		Histogram /b=1 ROI, Histo	// Histogram it
		wavestats /q/m=1 histo				// compute the wavestats... includes V_maxLoc -> the x location of the largest bin
		//if(V_MaxLoc>0)
		F0 = V_MaxLoc
		//	print F0
		//else
		//	F0 = 1
		//endif
		DFF0_pop[][i]=(ROI[p] -F0)/ F0
	endfor
	
	Duplicate/O DFF0_pop $DFFOwavename
	Duplicate/O regMovie $squarewavename
	Duplicate/O sqROIwave $ROIMaskwavename

	KillWaves/Z/F SurroundnPixels, ROI, ROInPixels, sqROIWave, ExclusionMask, regMovie, BKGCorr_pop,ROI_pop, BKG_pop
	KillWaves/Z/F ROInPixels, ROIcentres, ROIcentres, surround_ROI, DFF0_pop

End
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

Function vonMises(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = A*exp(k*(cos(2*(x - T)) - 1))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = k
	//CurveFitDialog/ w[2] = T

	return w[0]*exp(w[1]*(cos(2*(x - w[2])) - 1))
End