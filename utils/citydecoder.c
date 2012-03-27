/****************************************************************************
*                CITYDECODER: A SPLAT! File Conversion Utility              *
*                  Copyright John A. Magliacane, KD2BD 2002                 *
*                         Last update: 08-Feb-2009                          *
*****************************************************************************
*                                                                           *
* This utility reads ASCII Metadata Cartographic Boundary Files available   *
* through the U.S. Census Bureau, and generates a lists of cities, states,  *
* and counties along with the latitude and longitude corresponding to their *
* geographic centers.  Such data may be (optionally) sorted and written to  *
* files for use with SPLAT! software.  This utility takes as an argument,   *
* two-letter prefix plus the FIPS code for the state being processed (ie:   *
* "citydecoder pl34" will read files "pl34_d00.dat" and "pl34_d00a.dat",    *
* and produce a list of city names and geographical coordinates for the     *
* state of New Jersey.                                                      *
*                                                                           *
*****************************************************************************
*                                                                           *
* This program is free software; you can redistribute it and/or modify it   *
* under the terms of the GNU General Public License as published by the     *
* Free Software Foundation; either version 2 of the License or any later    *
* version.                                                                  *
*                                                                           *
* This program is distributed in the hope that it will useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     *
* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License     *
* for more details.                                                         *
*                                                                           *
*****************************************************************************
*          To compile: gcc -Wall -O6 -s citydecoder.c -o citydecoder        *
*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(argc,argv)
char argc, *argv[];
{
	int x, y, z, n;
	long attributefile_id, coordinatefile_id;
	char string[80], name[80], attributefilename[15], coordinatefilename[15], *s=NULL;
	double lat, lon;
	FILE *attributefile=NULL, *coordinatefile=NULL;

	if (argc==1)
	{
		fprintf(stderr,"\n*** Usage: citydecoder pl34 pl36 pl42 | sort > outputfile\n\n");
		exit(1);
	}

	for (z=1; z<argc; z++)
	{
		sprintf(attributefilename,"%s_d00a.dat",argv[z]);
		sprintf(coordinatefilename,"%s_d00.dat",argv[z]);

		attributefile=fopen(attributefilename,"r");
		coordinatefile=fopen(coordinatefilename,"r");

		if (attributefile!=NULL && coordinatefile!=NULL)
		{
			/* Skip First ASCII File Record (ID=0) */

			for (x=0; x<7; x++)
				s=fgets(string,80,attributefile);

			do
			{
				string[0]=0;
				n=fscanf(coordinatefile,"%ld", &coordinatefile_id);

				if (coordinatefile_id!=-99999)
				{
					name[0]=0;

					n=fscanf(coordinatefile,"%lf %lf",&lon, &lat);

					/* Read ID Number From Attribute File */

					s=fgets(string,80,attributefile);
					n=sscanf(string,"%ld",&attributefile_id);


					/* Skip Two Strings in Attribute File */

					s=fgets(string,80,attributefile);
					s=fgets(string,80,attributefile);


					/* Read City Name From Attribute File */

					s=fgets(string,80,attributefile);

					/* Strip "quote" characters from name */

					for (x=2, y=0; string[x]!='"' && string[x]!=0; x++, y++)
						name[y]=string[x];

					name[y]=0;

					/* Skip Two Strings in Attribute File */

					s=fgets(string,80,attributefile);
					s=fgets(string,80,attributefile);

					/* Skip blank line between records */

					s=fgets(string,80,attributefile);

					if (name[0]!=0 && name[0]!=' ' && feof(attributefile)==0 && attributefile_id==coordinatefile_id)
						printf("%s, %f, %f\n",name,lat,-lon);
				}


				/* Read to the end of the current coordinatefile record */

				do
				{
					string[0]=0;
					n=fscanf(coordinatefile,"%s",string);

				} while (strncmp(string,"END",3)!=0 && feof(coordinatefile)==0);

			} while (feof(coordinatefile)==0);

			fclose(attributefile);
			fclose(coordinatefile);
		}

		else
		{
			/* Houston, We Have A Problem... */

			fprintf(stderr,"%c\n",7);

			if (coordinatefile==NULL)
				fprintf(stderr,"*** Error opening coordinate file: \"%s\"!\n",coordinatefilename);

			if (attributefile==NULL)
				fprintf(stderr,"*** Error opening attribute file : \"%s\"!\n",attributefilename);
			fprintf(stderr,"\n");
		}
	}
		
	exit(0);
}
