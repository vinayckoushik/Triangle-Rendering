#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* 3x3 * 3x1 Matrix Multiplication
void matrix_mult(float a[3],float b[3][3],float c[3])
{

	c[0]=a[0]*b[0][0]+a[1]*b[1][0]+a[2]*b[2][0];
	c[1]=a[0]*b[0][1]+a[1]*b[1][1]+a[2]*b[2][1];
	c[2]=a[0]*b[0][2]+a[1]*b[1][2]+a[2]*b[2][2];

}

/*Cross Product*/
void cross_prod(float ax,float ay,float az,float bx,float by,float bz,float c[3])
{
	c[0] = (ay*bz)-(az*by);
	c[1] = (az*bx)-(ax*bz);
	c[2] = (ax*by)-(ay*bx);

}

/*Dot Product*/
void dot_prod(float ax,float ay,float az,float bx,float by,float bz,float *c)
{
	*c=(ax*bx)+(ay*by)+(az*bz);
}

int main(int argc,char *argv[])
{
	FILE *fpt,*fpo;
	int i=0,j=0,k=0;
	fpt=fopen(argv[1],"rb");
	fpo=fopen("Triangle.ppm","wb");
	unsigned char C[100];
	unsigned char al;
	int total,face,xx;

	/*PARSING HEADER*/

	while(strcmp(C,"end_header")!=0)
	{
		if(strcmp(C,"vertex")==0)
			fscanf(fpt,"%d",&total);
		else if(strcmp(C,"face")==0)
			fscanf(fpt,"%d",&face);
		fscanf(fpt,"%s",C);
		i++;
	}

	/*DECLARATION*/

	float x[total],y[total],z[total],temp=0,min_x=0,min_y=0,min_z=0,max_x=0,max_y=0,max_z=0;
	float center_x,center_y,center_z;
	float extent_x,extent_y,extent_z,extent;
	float camera_x=1,camera_y=0,camera_z=0;
	float up_x=0,up_y=0,up_z=1;
	float rx[3][3],ry[3][3],rz[3][3];
	float theta_x,theta_y,theta_z;
	float m1[3],m2[3];
	float left[3],right[3],top[3],bottom[3],topleft[3];
	float a;
	int r,c,COLS=256,ROWS=256;
	float image[ROWS*COLS][3];
	unsigned char ip[ROWS*COLS];
	float A,B,C1,D,eqn[3];
	float n,d;
	float intersect[3];
	float dot1,dot2,dot3;
	float temp1[3],temp2[3];
	float camera[3];
	float up[3];

	int xf[face],yf[face],zf[face];

	fprintf(fpo,"P5 256 256 255\n");

	/*READING VERTICES AND FACES*/

	for(i=0;i<total;i++)
	{
		fscanf(fpt,"%f %f %f",&x[i],&y[i],&z[i]);
	}
	
	for(i=0;i<face;i++)
	{
		fscanf(fpt,"%d %d %d %d",&xx,&xf[i],&yf[i],&zf[i]);
	}
	
	min_x=x[0];
	min_y=y[0];
	min_z=z[0];

	/*FINDING MAX AND MIN*/

	for(i=0;i<total;i++)
	{
		if(x[i]<min_x)
			min_x=x[i];
		else if(y[i]<min_y)
			min_y=y[i];
		else if(z[i]<min_z)
			min_z=z[i];
		else if(x[i]>max_x)
			max_x=x[i];
		else if(y[i]>max_y)
			max_y=y[i];
		else if(z[i]>max_z)
			max_z=z[i];
	}

	center_x=(min_x+max_x)/2;
	center_y=(min_y+max_y)/2;
	center_z=(min_z+max_z)/2;

	extent_x=max_x-min_x;
	extent_y=max_y-min_y;
	extent_z=max_z-min_z;

	if(extent_x > extent_y && extent_x > extent_z)
		extent=extent_x;
	else if(extent_y > extent_x && extent_y > extent_z)
		extent=extent_y;
	else if(extent_z > extent_y && extent_z > extent_x)
		extent=extent_z;

	/*READING THETA*/

	theta_x=atoi(argv[2]);
	theta_y=atoi(argv[3]);
	theta_z=atoi(argv[4]);

	/*DEFINING ROTATION MATRICES*/

	rx[0][0]=1;
	rx[0][1]=0;
	rx[0][2]=0;
	rx[1][0]=0;
	rx[1][1]=cos(M_PI*theta_x/180.0);
	rx[1][2]=-sin(M_PI*theta_x/180.0);
	rx[2][0]=0;
	rx[2][1]=sin(M_PI*theta_x/180.0);
	rx[2][2]=cos(M_PI*theta_x/180.0);

	ry[0][0]=cos(M_PI*theta_y/180.0);
	ry[0][1]=0;
	ry[0][2]=sin(M_PI*theta_y/180.0);
	ry[1][0]=0;
	ry[1][1]=1;
	ry[1][2]=0;
	ry[2][0]=-sin(M_PI*theta_y/180.0);
	ry[2][1]=0;
	ry[2][2]=cos(M_PI*theta_y/180.0);

	rz[0][0]=cos(M_PI*theta_z/180.0);
	rz[0][1]=-sin(M_PI*theta_z/180.0);
	rz[0][2]=0;
	rz[1][0]=sin(M_PI*theta_z/180.0);
	rz[1][1]=cos(M_PI*theta_z/180.0);
	rz[1][2]=0;
	rz[2][0]=0;
	rz[2][1]=0;
	rz[2][2]=1;

	camera[0]=camera_x;
	camera[1]=camera_y;
	camera[2]=camera_z;

	/*MATRIX MULTIPLICATION OF CAMERA VECTOR*/

	matrix_mult(camera,rx,m1);
	matrix_mult(m1,ry,m2);
	matrix_mult(m2,rz,camera);

	up[0]=up_x;
	up[1]=up_y;
	up[2]=up_z;

	/*ROTATING UP VECTOR*/

	matrix_mult(up,rx,m1);
	matrix_mult(m1,ry,m2);
	matrix_mult(m2,rz,up);

	camera_x=camera[0];
	camera_y=camera[1];
	camera_z=camera[2];

	up_x=up[0];
	up_y=up[1];
	up_z=up[2];

	/*SCALING AND MULTIPLYING*/

	camera_x=1.5*extent*camera_x +center_x;
	camera_y=1.5*extent*camera_y +center_y;
	camera_z=1.5*extent*camera_z +center_z;

	/*FINDING LEFT*/

	cross_prod(up_x,up_y,up_z,(center_x-camera_x),(center_y-camera_y),(center_z-camera_z),left);

	a=sqrt(pow(left[0],2)+pow(left[1],2)+pow(left[2],2));

	left[0]=(extent/(2*a))*left[0]+center_x;
	left[1]=(extent/(2*a))*left[1]+center_y;
	left[2]=(extent/(2*a))*left[2]+center_z;

	cross_prod((center_x-camera_x),(center_y-camera_y),(center_z-camera_z),up_x,up_y,up_z,right);
	
	right[0]=(extent/(2*a))*right[0]+center_x;
	right[1]=(extent/(2*a))*right[1]+center_y;
	right[2]=(extent/(2*a))*right[2]+center_z;

	top[0]=(extent/(2.0))*up_x+center_x;
	top[1]=(extent/(2.0))*up_y+center_y;
	top[2]=(extent/(2.0))*up_z+center_z;

	bottom[0]=(-extent/(2.0))*up_x+center_x;
	bottom[1]=(-extent/(2.0))*up_y+center_y;
	bottom[2]=(-extent/(2.0))*up_z+center_z;

	topleft[0]=(extent/(2.0))*up_x+left[0];
	topleft[1]=(extent/(2.0))*up_y+left[1];
	topleft[2]=(extent/(2.0))*up_z+left[2];

	for(r=0;r<(ROWS);r++)
	{
		for(c=0;c<(COLS);c++)
		{
			image[r*COLS+c][0]=topleft[0]+((float) c/(float)(COLS-1))*(right[0]-left[0])+((float) r/((float)(ROWS-1)))*(bottom[0]-top[0]);
			image[r*COLS+c][1]=topleft[1]+((float) c/(float)(COLS-1))*(right[1]-left[1])+((float) r/((float)(ROWS-1)))*(bottom[1]-top[1]);
			image[r*COLS+c][2]=topleft[2]+((float) c/(float)(COLS-1))*(right[2]-left[2])+((float) r/((float)(ROWS-1)))*(bottom[2]-top[2]);
		}
	}

	/*CLEARING O/P IMAGE*/

	for(i=0;i<ROWS*COLS;i++)
		ip[i]=0;
	
	/*FILLING PIXELS FOR EACH TRIANGLE*/

	for(k=0;k<face;k++)
	{

		cross_prod((x[yf[k]]-x[xf[k]]),(y[yf[k]]-y[xf[k]]),(z[yf[k]]-z[xf[k]]),(x[zf[k]]-x[xf[k]]),(y[zf[k]]-y[xf[k]]),(z[zf[k]]-z[xf[k]]),eqn);
		A=eqn[0];
		B=eqn[1];
		C1=eqn[2];

		dot_prod(-A,-B,-C1,x[xf[k]],y[xf[k]],z[xf[k]],&D);
		
		dot_prod(-A,-B,-C1,(camera_x),(camera_y),(camera_z),&n);

		n=n-D;


		for(i=0;i<ROWS*COLS;i++)
		{
			if(i%65535==0)
				printf("%d ",j++);

			dot_prod(A,B,C1,(image[i][0]-camera_x),(image[i][1]-camera_y),(image[i][2]-camera_z),&d);

			if(fabs(d)<0.0001)
				continue;

			intersect[0]=camera_x+(n/d)*(image[i][0]-camera_x);
			intersect[1]=camera_y+(n/d)*(image[i][1]-camera_y);
			intersect[2]=camera_z+(n/d)*(image[i][2]-camera_z);

			cross_prod((x[zf[k]]-x[xf[k]]),(y[zf[k]]-y[xf[k]]),(z[zf[k]]-z[xf[k]]),(x[yf[k]]-x[xf[k]]),(y[yf[k]]-y[xf[k]]),(z[yf[k]]-z[xf[k]]),temp1);
			cross_prod((intersect[0]-x[xf[k]]),(intersect[1]-y[xf[k]]),(intersect[2]-z[xf[k]]),(x[yf[k]]-x[xf[k]]),(y[yf[k]]-y[xf[k]]),(z[yf[k]]-z[xf[k]]),temp2);
			dot_prod(temp1[0],temp1[1],temp1[2],temp2[0],temp2[1],temp2[2],&dot1);

			cross_prod((x[xf[k]]-x[yf[k]]),(y[xf[k]]-y[yf[k]]),(z[xf[k]]-z[yf[k]]),(x[zf[k]]-x[yf[k]]),(y[zf[k]]-y[yf[k]]),(z[zf[k]]-z[yf[k]]),temp1);
			cross_prod((intersect[0]-x[yf[k]]),(intersect[1]-y[yf[k]]),(intersect[2]-z[yf[k]]),(x[zf[k]]-x[yf[k]]),(y[zf[k]]-y[yf[k]]),(z[zf[k]]-z[yf[k]]),temp2);
			dot_prod(temp1[0],temp1[1],temp1[2],temp2[0],temp2[1],temp2[2],&dot2);

			cross_prod((x[yf[k]]-x[zf[k]]),(y[yf[k]]-y[zf[k]]),(z[yf[k]]-z[zf[k]]),(x[xf[k]]-x[zf[k]]),(y[xf[k]]-y[zf[k]]),(z[xf[k]]-z[zf[k]]),temp1);
			cross_prod((intersect[0]-x[zf[k]]),(intersect[1]-y[zf[k]]),(intersect[2]-z[zf[k]]),(x[xf[k]]-x[zf[k]]),(y[xf[k]]-y[zf[k]]),(z[xf[k]]-z[zf[k]]),temp2);
			dot_prod(temp1[0],temp1[1],temp1[2],temp2[0],temp2[1],temp2[2],&dot3);

			if((dot1>0) && (dot2>0) && (dot3>0))
			{
				ip[i]=155+(k%100);
			}
		}
		
	}
	
	fwrite(ip,1,ROWS*COLS,fpo);
	

	fclose(fpt);
	fclose(fpo);
}