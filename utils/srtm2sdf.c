/**************************************************************\
 **     Created originally by Jonathan Naylor, G4KLX.        **
 **     Later embellished by John Magliacane, KD2BD to       **
 **     detect and handle voids found in the SRTM data,      **
 **     SRTM-3 data in .BIL and .HGT format, and high        **
 **     resolution SRTM-1 one arc-second topography data.    **
 **************************************************************
 **                    Compile like this:                    **
 **      cc -Wall -O3 -s -lbz2 srtm2sdf.c -o srtm2sdf        **
 **              Last modification: 12-Jan-2019             **
\**************************************************************/ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <bzlib.h>

#define BZBUFFER 65536

char	sdf_filename[30], sdf_path[255], replacement_flag, opened=0,
	hgt=0, bil=0;
char colons=0;

int	srtm[3601][3601], usgs[1201][1201], max_north, max_west, n,
	min_north, min_west, merge=0, min_elevation, bzerror, ippd, mpi;

int ReadSRTM(char *filename)
{
	int x, y, byte=0;
	unsigned char error, buffer[2];
	char north[3], west[4], *base=NULL, blw_filename[255];
	double cell_size, deg_north, deg_west;
	FILE *fd=NULL;

	if (strstr(filename, ".zip")!=NULL)
	{
		fprintf(stderr, "*** Error: \"%s\" must be uncompressed\n",filename);
		return -1;

	}

	if (strstr(filename, ".tgz")!=NULL)
	{
		fprintf(stderr, "*** Error: \"%s\" must be uncompressed\n",filename);
		return -1;

	}

	if ((strstr(filename, ".hgt")==NULL) && (strstr(filename, ".bil")==NULL))
	{
		fprintf(stderr, "*** Error: \"%s\" does not have the correct extension (.hgt or .bil)\n",filename);
		return -1;
	}

	if (strstr(filename, ".hgt")!=NULL)
		hgt=1;

	if (strstr(filename, ".bil")!=NULL)
		bil=1;

	base=strrchr(filename, '/');

	if (base==NULL)
		base=filename;
	else
		base+=1;

	if (hgt)
	{
		/* We obtain coordinates from the base of the .HGT filename */

		north[0]=base[1];
		north[1]=base[2];
		north[2]=0;

		west[0]=base[4];
		west[1]=base[5];
		west[2]=base[6];
		west[3]=0;

		if ((base[0]!='N' && base[0]!='S') || (base[3]!='W' && base[3]!='E'))
		{
			fprintf(stderr, "*** Error: \"%s\" doesn't look like a valid .hgt SRTM filename.\n", filename);
			return -1;
		}

		max_west=atoi(west);

		if (base[3]=='E')
			max_west=360-max_west;

		min_west=max_west-1;

		if (max_west==360)
			max_west=0;

		if (base[0]=='N')
			min_north=atoi(north);
		else
			min_north=-atoi(north);

		max_north=min_north+1;
	}

	if (bil)
	{
		/* We obtain .BIL file coordinates
		   from the corresponding .BLW file */

		strncpy(blw_filename,filename,250);
		x=(int)strlen(filename);

		if (x>3)
		{
			blw_filename[x-2]='l';
			blw_filename[x-1]='w';
			blw_filename[x]=0;

			fd=fopen(blw_filename,"rb");

			if (fd!=NULL)
			{
				fscanf(fd,"%lf",&cell_size);

				if ((cell_size<0.0008) || (cell_size>0.0009))
				{
					printf("\n*** .BIL file's cell size is incompatible with SPLAT!!\n");
					exit(1);
				}

				fscanf(fd,"%lf",&deg_west);
				fscanf(fd,"%lf",&deg_west);
				fscanf(fd,"%lf",&deg_west);

				fscanf(fd,"%lf",&deg_west);

				fscanf(fd,"%lf",&deg_north);

				fclose(fd);
			}

			min_north=(int)(deg_north);
			max_north=max_north+1;

			if (deg_west<0.0)
				deg_west=-deg_west;
			else
				deg_west=360.0-deg_west;

			min_west=(int)(deg_west);

			if (min_west==360)
				min_west=0;

			max_west=min_west+1;
		}
	}

	fd=fopen(filename, "rb");

	if (!fd)
	{
		fprintf(stderr, "*** Error: Cannot open \"%s\"\n", filename);
		return -1;
	}

	fread(&buffer,1,2,fd);

	if ((buffer[0]=='P') && (buffer[1]=='K'))
	{
		fprintf(stderr, "*** Error: \"%s\" still appears to be compressed!\n",filename);
		fclose(fd);
		return -1;
	}

	fseek(fd,0L,SEEK_SET);

    if (colons) {
		sprintf(sdf_filename, "%d:%d:%d:%d%s.sdf", min_north, max_north, min_west, max_west, (ippd==3600?"-hd":""));
    } else {
		sprintf(sdf_filename, "%d_%d_%d_%d%s.sdf", min_north, max_north, min_west, max_west, (ippd==3600?"-hd":""));
    }

	error=0;
	replacement_flag=0;

	printf("Reading %s... ", filename);
	fflush(stdout);

	for (x=0; (x<=ippd && error==0); x++)
		for (y=0; (y<=ippd && error==0); y++)
		{

			if (fread(&buffer, 1, 2, fd)==2)
			{
				if (bil)
				{
					/* "little-endian" structure */

					byte=buffer[0]+(buffer[1]<<8);

					if (buffer[1]&128)
						byte-=0x10000;
				}

				if (hgt)
				{
					/* "big-endian" structure */

					byte=buffer[1]+(buffer[0]<<8);

					if (buffer[0]&128)
						byte-=0x10000;
				}

				/* Flag problem elevations here */

				if (byte<-32768)
					byte=-32768;

				if (byte>32767)
					byte=32767;

				if (byte<=min_elevation)
					replacement_flag=1;

				srtm[x][y]=byte;
			}

			else
				error=1;
		}

	if (error)
	{
		fprintf(stderr,"\n*** Error: Premature EOF detected while reading \"%s\"!  :-(\n",filename);
	}

	fclose(fd);

	return 0;
}

int LoadSDF_SDF(char *name)
{
	/* This function reads uncompressed
	   SPLAT Data Files (.sdf) into memory. */

	int x, y, dummy;
	char sdf_file[255], path_plus_name[512];
	FILE *infile;

	for (x=0; name[x]!='.' && name[x]!=0 && x<250; x++)
		sdf_file[x]=name[x];

	sdf_file[x]='.';
	sdf_file[x+1]='s';
	sdf_file[x+2]='d';
	sdf_file[x+3]='f';
	sdf_file[x+4]=0;

	strncpy(path_plus_name,sdf_path,255);
	strncat(path_plus_name,sdf_file,254);

	infile=fopen(path_plus_name,"rb");

	if (infile==NULL)
		return 0;

	fscanf(infile,"%d", &dummy);
	fscanf(infile,"%d", &dummy);
	fscanf(infile,"%d", &dummy);
	fscanf(infile,"%d", &dummy);

	printf("\nReading %s... ",path_plus_name);
	fflush(stdout);

	for (x=0; x<1200; x++)
		for (y=0; y<1200; y++)
			fscanf(infile,"%d",&usgs[x][y]);

	fclose(infile);

	return 1;
}

char *BZfgets(BZFILE *bzfd, unsigned length)
{
	/* This function returns at most one less than 'length' number
	   of characters from a bz2 compressed file whose file descriptor
	   is pointed to by *bzfd.  In operation, a buffer is filled with
	   uncompressed data (size = BZBUFFER), which is then parsed
	   and doled out as NULL terminated character strings every time
	   this function is invoked.  A NULL string indicates an EOF
	   or error condition. */

	static int x, y, nBuf;
	static char buffer[BZBUFFER+1], output[BZBUFFER+1];
	char done=0;

	if (opened!=1 && bzerror==BZ_OK)
	{
		/* First time through.  Initialize everything! */

		x=0;
		y=0;
		nBuf=0;
		opened=1;
		output[0]=0;
	}

	do
	{
		if (x==nBuf && bzerror!=BZ_STREAM_END && bzerror==BZ_OK && opened)
		{
			/* Uncompress data into a static buffer */

			nBuf=BZ2_bzRead(&bzerror, bzfd, buffer, BZBUFFER);
			buffer[nBuf]=0;
			x=0;
		}

		/* Build a string from buffer contents */

		output[y]=buffer[x];

		if (output[y]=='\n' || output[y]==0 || y==(int)length-1)
		{
			output[y+1]=0;
			done=1;
			y=0;
		}

		else
			y++;
		x++;

	} while (done==0);

	if (output[0]==0)
		opened=0;

	return (output);
}

int LoadSDF_BZ(char *name)
{
	/* This function reads .bz2 compressed
	   SPLAT Data Files into memory. */

	int x, y, dummy;
	char sdf_file[255], path_plus_name[255];
	FILE *fd;
	BZFILE *bzfd;

	for (x=0; name[x]!='.' && name[x]!=0 && x<247; x++)
		sdf_file[x]=name[x];

	sdf_file[x]='.';
	sdf_file[x+1]='s';
	sdf_file[x+2]='d';
	sdf_file[x+3]='f';
	sdf_file[x+4]='.';
	sdf_file[x+5]='b';
	sdf_file[x+6]='z';
	sdf_file[x+7]='2';
	sdf_file[x+8]=0;

	strncpy(path_plus_name,sdf_path,255);
	strncat(path_plus_name,sdf_file,254);

	fd=fopen(path_plus_name,"rb");
	bzfd=BZ2_bzReadOpen(&bzerror,fd,0,0,NULL,0);

	if (fd!=NULL && bzerror==BZ_OK)
	{
		printf("\nReading %s... ",path_plus_name);
		fflush(stdout);

		sscanf(BZfgets(bzfd,255),"%d",&dummy);
		sscanf(BZfgets(bzfd,255),"%d",&dummy);
		sscanf(BZfgets(bzfd,255),"%d",&dummy);
		sscanf(BZfgets(bzfd,255),"%d",&dummy);
	
		for (x=0; x<1200; x++)
			for (y=0; y<1200; y++)
				sscanf(BZfgets(bzfd,20),"%d",&usgs[x][y]);

		fclose(fd);

		BZ2_bzReadClose(&bzerror,bzfd);

		return 1;
	}

	else
		return 0;
}

char LoadSDF(char *name)
{
	/* This function loads the requested SDF file from the filesystem.
	   First, it tries to invoke the LoadSDF_SDF() function to load an
	   uncompressed SDF file (since uncompressed files load slightly
	   faster).  Failing that, it tries to load a compressed SDF file
	   by invoking the LoadSDF_BZ() function. */

	int return_value=-1;

	/* Try to load an uncompressed SDF first. */

	return_value=LoadSDF_SDF(name);

	/* If that fails, try loading a compressed SDF. */

	if (return_value==0 || return_value==-1)
		return_value=LoadSDF_BZ(name);

	return return_value;
}

int ReadUSGS()
{
	char usgs_filename[15];

	/* usgs_filename is a minimal filename ("40:41:74:75").
	   Full path and extentions are added later though
	   subsequent function calls. */

	sprintf(usgs_filename, "%d:%d:%d:%d", min_north, max_north, min_west, max_west);

	return (LoadSDF(usgs_filename));
}

void average_terrain(y,x,z)
int x, y, z;
{
	long accum;
	int temp=0, count, bad_value;
	double average;

	bad_value=srtm[y][x];

	accum=0L;
	count=0;

	if (y>=2)
	{
		temp=srtm[y-1][x];

		if (temp>bad_value)
		{
			accum+=temp;
			count++;
		}
	}

	if (y<=mpi)
	{
		temp=srtm[y+1][x];

		if (temp>bad_value)
		{
			accum+=temp;
			count++;
		}
	}

	if ((y>=2) && (x<=(mpi-1)))
	{
		temp=srtm[y-1][x+1];

		if (temp>bad_value)
		{
			accum+=temp;
			count++;
		}
	}

	if (x<=(mpi-1))
	{
		temp=srtm[y][x+1];

		if (temp>bad_value)
		{
			accum+=temp;
			count++;
		}
	}

	if ((x<=(mpi-1)) && (y<=mpi))
	{
		temp=srtm[y+1][x+1];

		if (temp>bad_value)
		{
			accum+=temp;
			count++;
		}
	}

	if ((x>=1) && (y>=2)) 
	{
		temp=srtm[y-1][x-1];

		if (temp>bad_value)
		{
			accum+=temp;
			count++;
		}
	}

	if (x>=1)
	{
		temp=srtm[y][x-1];

		if (temp>bad_value)
		{
			accum+=temp;
			count++;
		}
	}

	if ((y<=mpi) && (x>=1))
	{
		temp=srtm[y+1][x-1];

		if (temp>bad_value)
		{
			accum+=temp;
			count++;
		}
	}

	if (count!=0)
	{
		average=(((double)accum)/((double)count));
		temp=(int)(average+0.5);
	}

	if (temp>min_elevation)
		srtm[y][x]=temp;
	else
		srtm[y][x]=min_elevation;
}

void WriteSDF(char *filename)
{
	/* Like the HGT files, the extreme southwest corner
	 * provides the point of reference for the SDF file.
	 * The SDF file extends from min_north degrees to
	 * the south to min_north+(mpi/ippd) degrees to
	 * the north, and max_west degrees to the west to
	 * max_west-(mpi/ippd) degrees to the east.  The
	 * overlapping edge redundancy present in the HGT
	 * and earlier USGS files is not necessary, nor
	 * is it present in SDF files.
	 */

	int x, y, byte, last_good_byte=0;
	FILE *outfile;

	printf("\nWriting %s... ", filename);
	fflush(stdout);

	outfile=fopen(filename,"wb");

	fprintf(outfile, "%d\n%d\n%d\n%d\n", max_west, min_north, min_west, max_north);

	for (y=ippd; y>=1; y--)		/* Omit the northern most edge */
		for (x=mpi; x>=0; x--) /* Omit the eastern most edge  */
		{
			byte=srtm[y][x];

			if (byte>min_elevation)
				last_good_byte=byte;

			if (byte<min_elevation)
			{
				if (merge)
				{
					if (ippd==3600)
						fprintf(outfile,"%d\n",usgs[1200-(y/3)][1199-(x/3)]);
					else
						fprintf(outfile,"%d\n",usgs[1200-y][1199-x]);
				}

				else
				{
					average_terrain(y,x,last_good_byte);
					fprintf(outfile,"%d\n",srtm[y][x]);
				}
			}
			else
				fprintf(outfile,"%d\n",byte);
		}

	printf("Done!\n");

	fclose(outfile);
}

int main(int argc, char **argv)
{
	int x, y, z=0;
	char *env=NULL;

	if (argc==1 || (argc==2 && strncmp(argv[1],"-h",2)==0))
	{
        fprintf(stderr, "\nsrtm2sdf: Generates SPLAT! elevation data files from unzipped\nSRTM-1 or SRTM-3 elevation data files, and replaces SRTM data voids with\nelevation data from older usgs2sdf derived SDF files.\n\n");

		fprintf(stderr, "\tAvailable Options...\n\n");
		fprintf(stderr, "\t-c Use colons in the filename instead of x's\n\t    This mimics older behavior, but generates files that cannot be used on Windows filesystems.\n\n");
		fprintf(stderr, "\t-d directory path of usgs2sdf derived SDF files\n\t    (overrides path in ~/.splat_path file)\n\n");
		fprintf(stderr, "\t-hd Invoke in High Definition mode, working on SRTM-1 files\n\t    Without this, it expects the input files to be SRTM-3 files.\n\t    You can also invoke the binary as %s-hd\n\n", argv[0]);
		fprintf(stderr, "\t-n elevation limit (in meters) below which SRTM data is\n\t    replaced by USGS-derived .sdf data (default = 0 meters)\n\n");
		fprintf(stderr, "Examples: %s N40W074.hgt\n",argv[0]);
		fprintf(stderr, "          %s -d /cdrom/sdf N40W074.hgt\n",argv[0]);
		fprintf(stderr, "          %s -d /dev/null N40W074.hgt  (prevents data replacement)\n",argv[0]);
		fprintf(stderr, "          %s -n -5 N40W074.hgt\n\n",argv[0]);

		return 1;
	}

	if (strstr(argv[0], "-hd")!=NULL) {
		ippd=3600;	/* High Definition (1 arc-sec) Mode */
	} else {
		ippd=1200;	/* Standard Definition (3 arc-sec) Mode */
	}

	y=argc-1;

	min_elevation=0;

	for (x=1, z=0; x<=y; x++)
	{
		if (strcmp(argv[x],"-c")==0)
		{
			z=x;
            colons=1;
		}

		if (strcmp(argv[x],"-d")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
				strncpy(sdf_path,argv[z],253);
		}

		if (strcmp(argv[x],"-hd")==0)
		{
			z=x;
            ippd=3600;	/* High Definition (1 arc-sec) Mode */
		}

		if (strcmp(argv[x],"-n")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0])
			{
				sscanf(argv[z],"%d",&min_elevation);

				if (min_elevation<-32767)
					min_elevation=0;
			}			 
		}
	}

	mpi=ippd-1;		/* Maximum pixel index per degree */

	/* If no SDF path was specified on the command line (-d), check
	   for a path specified in the $HOME/.splat_path file.  If the
	   file is not found, then sdf_path[] remains NULL, and a data
	   merge will not be attempted if voids are found in the SRTM file. */

#ifndef _WIN32
	if (sdf_path[0]==0)
	{
		env=getenv("HOME");

		sprintf(string,"%s/.splat_path",env);

		fd=fopen(string,"r");

		if (fd!=NULL)
		{
			fgets(string,253,fd);

			/* Remove <CR> and/or <LF> from string */

			for (x=0; string[x]!=13 && string[x]!=10 && string[x]!=0 && x<253; x++);
			string[x]=0;

			strncpy(sdf_path,string,253);

			fclose(fd);
		}
	}
#endif

	/* Ensure a trailing '/' is present in sdf_path */

	if (sdf_path[0])
	{
		x=(int)strlen(sdf_path);

		if (sdf_path[x-1]!='/' && x!=0)
		{
			sdf_path[x]='/';
			sdf_path[x+1]=0;
		}
	}

	if (ReadSRTM(argv[z+1])==0)
	{
		if (replacement_flag && sdf_path[0])
			merge=ReadUSGS();

		WriteSDF(sdf_filename);
	}

	return 0;
}

