	//the function to calculate free engergy at temperature T(dg),with dg 37(data) and dh 37(dhdata) 

inline short Tscale(float T,short dG, short dH)
{
	//return (short)floor((float)(dH)-((float)(dH-dG)*T/310.15)+0.5);
//	if (dG==dH)
//	return dG;
//	else
	return (short)floor((float)dG*T/310.15 -(float)dH*T*(T-310.15)/(310.15*310.15)+0.5);
}

void dG_T (float T, datatable &data, datatable &dhdata, datatable &dg)
{
	int a,b,c,d,e,f,g,h;
	
	//if the temperature is 37 degree dg equals data
	dg=data;
	//calculated with T
	dg.prelog=data.prelog*T/310.15;
    dg.RT=data.RT*T/310.15;
	 for(a=1;a<=data.numoftloops;a++)
   	 (dg.tloop[a][1])=Tscale(T,data.tloop[a][1],dhdata.tloop[a][1]);
	 for(a=1;a<=data.numoftriloops;a++)
	 (dg.triloop[a][1])=Tscale(T,data.triloop[a][1],dhdata.triloop[a][1]);
	 for(a=1;a<=data.numofhexaloops;a++)
 	 (dg.hexaloop[a][1])=Tscale(T,data.hexaloop[a][1],dhdata.hexaloop[a][1]);		
  			    
	
	for(a=1;a<5;a++)
			dg.poppen[a]=Tscale(T,data.poppen[a],dhdata.poppen[a]);
		
			dg.maxpen=Tscale(T,data.maxpen,dhdata.maxpen);
	for(a=1;a<7;a++)
		        dg.eparam[a]=Tscale(T,data.eparam[a],dhdata.eparam[a]);
		        dg.eparam[10]=Tscale(T,data.eparam[10],dhdata.eparam[10]);

// apply internal initiation for internal loop larger than 5
	for (a=6;a<=30;a++) {
		dhdata.inter[a]=ARGV[4];
	}

	for(a=1;a<31;a++)

	{		 
		dg.inter[a]=Tscale(T,data.inter[a],dhdata.inter[a]);
		dg.bulge[a]=Tscale(T,data.bulge[a],dhdata.bulge[a]);
		dg.hairpin[a]=Tscale(T,data.hairpin[a],dhdata.hairpin[a]);
	

	}

	dg.auend=Tscale(T,data.auend,dhdata.auend);
	dg.gubonus=Tscale(T,data.gubonus,dhdata.gubonus);
	dg.cint=Tscale(T,data.cint,dhdata.cint);
	dg.cslope=Tscale(T,data.cslope,dhdata.cslope);
	dg.c3=Tscale(T,data.c3,dhdata.c3);
	dg.efn2a=Tscale(T,data.efn2a,dhdata.efn2a);
	dg.efn2b=Tscale(T,data.efn2b,dhdata.efn2b);
	dg.efn2c=Tscale(T,data.efn2c,dhdata.efn2c);
	dg.init=Tscale(T,data.init,dhdata.init);
	dg.mlasym=Tscale(T,data.mlasym,dhdata.mlasym);
	dg.strain=Tscale(T,data.strain,dhdata.strain);
	dg.singlecbulge=Tscale(T,data.singlecbulge,dhdata.singlecbulge);


	for(a=0;a<6;a++){
		for(b=0;b<6;b++){
		for(c=0;c<6;c++){
		for(d=0;d<6;d++){
		for(e=0;e<6;e++){
		for(f=0;f<6;f++){
			(dg.iloop11[a][b][c][d][e][f])=Tscale(T,data.iloop11[a][b][c][d][e][f],dhdata.iloop11[a][b][c][d][e][f]);
		for(g=0;g<6;g++){
			(dg.iloop21[a][b][c][d][e][f][g])=Tscale(T,data.iloop21[a][b][c][d][e][f][g],dhdata.iloop21[a][b][c][d][e][f][g]);
			
		for(h=0;h<6;h++){
			(dg.iloop22[a][b][c][d][e][f][g][h])=Tscale(T,data.iloop22[a][b][c][d][e][f][g][h],dhdata.iloop22[a][b][c][d][e][f][g][h]);
			}
			}
			}
			}
		(dg.tstki1n[a][b][c][d])=Tscale(T,data.tstki1n[a][b][c][d],dhdata.tstki1n[a][b][c][d]);
		(dg.tstki23[a][b][c][d])=Tscale(T,data.tstki23[a][b][c][d],dhdata.tstki23[a][b][c][d]);	
		(dg.tstkm[a][b][c][d])=Tscale(T,data.tstkm[a][b][c][d],dhdata.tstkm[a][b][c][d]);
		(dg.tstkh[a][b][c][d])=Tscale(T,data.tstkh[a][b][c][d],dhdata.tstkh[a][b][c][d]);
		(dg.tstack[a][b][c][d])=Tscale(T,data.tstack[a][b][c][d],dhdata.tstack[a][b][c][d]);	
		(dg.coaxstack[a][b][c][d])=Tscale(T,data.coaxstack[a][b][c][d],dhdata.coaxstack[a][b][c][d]);
		(dg.tstackcoax[a][b][c][d])=Tscale(T,data.tstackcoax[a][b][c][d],dhdata.tstackcoax[a][b][c][d]);
		(dg.coax[a][b][c][d])=Tscale(T,data.coax[a][b][c][d],dhdata.coax[a][b][c][d]);	
		(dg.tstki[a][b][c][d])=Tscale(T,data.tstki[a][b][c][d],dhdata.tstki[a][b][c][d]);
		(dg.stack[a][b][c][d])=Tscale(T,data.stack[a][b][c][d],dhdata.stack[a][b][c][d]);	
		if(d==1||d==2)
		(dg.dangle[a][b][c][d])=Tscale(T,data.dangle[a][b][c][d],dhdata.dangle[a][b][c][d]);
			}
			}
			}
			}
			

}

