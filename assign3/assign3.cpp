/*
CSCI 480
Assignment 3 Raytracer

Name: <Your name here>
*/

#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <stdio.h>
#include <string>

#include "opencv2/core/core.hpp"
#include "opencv2/core/utility.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgcodecs/imgcodecs.hpp"

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);

#define vector(a,b,c) \
	(a)[0] = (b)[0] - (c)[0];	\
	(a)[1] = (b)[1] - (c)[1];	\
	(a)[2] = (b)[2] - (c)[2];

#define crossProduct(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

#define dotProduct(v,q) \
	((v)[0] * (q)[0] + \
	(v)[1] * (q)[1] + \
	(v)[2] * (q)[2])

#define scalarMultiply(a, b, c) \
	(a)[0] = (b)[0] * c; \
	(a)[1] = (b)[1] * c; \
	(a)[2] = (b)[2] * c;

// Triangle ray intersection algorithm from Moller and Trumbore
double rayIntersectsTriangle(float *p, float *d,
	Triangle triangle) {

	float e1[3], e2[3], h[3], s[3], q[3], v0[3], v1[3], v2[3];
	float a, f, u, v;

	v0[0] = triangle.v[0].position[0], v0[1] = triangle.v[0].position[1], v0[2] = triangle.v[0].position[2];
	v1[0] = triangle.v[1].position[0], v1[1] = triangle.v[1].position[1], v1[2] = triangle.v[1].position[2];
	v2[0] = triangle.v[2].position[0], v2[1] = triangle.v[2].position[1], v2[2] = triangle.v[2].position[2];

	vector(e1, v1, v0);
	vector(e2, v2, v0);

	crossProduct(h, d, e2);
	a = dotProduct(e1, h);

	if (a > -0.00001 && a < 0.00001)
		return(DBL_MAX);

	f = 1 / a;
	vector(s, p, v0);
	u = f * (dotProduct(s, h));

	if (u < 0.0 || u > 1.0)
		return(DBL_MAX);

	crossProduct(q, s, e1);
	v = f * dotProduct(d, q);

	if (v < 0.0 || u + v > 1.0)
		return(DBL_MAX);

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	double t = f * dotProduct(e2, q);

	if (t > 0.00001) // ray intersection
		return(t);

	else // this means that there is a line intersection
		 // but not a ray intersection
		return (DBL_MAX);

}

double rayIntersectsSphere(float *p, float *d, Sphere sphere) {
	double xc = sphere.position[0];
	double yc = sphere.position[1];
	double zc = sphere.position[2];
	double b = 2 * (d[0] * (p[0] - xc) + d[1] * (p[1] - yc) + d[2] * (p[2] - zc));
	double c = (p[0] - xc) * (p[0] - xc) + (p[1] - yc) * (p[1] - yc) + (p[2] - zc) * (p[2] - zc) - sphere.radius * sphere.radius;
	if (b * b - 4 * c >= 0) {
		double t0 = (-b + sqrt(b * b - 4 * c)) / 2;
		double t1 = (-b - sqrt(b * b - 4 * c)) / 2;
		double t = MIN(t0, t1);
		if (t > 0) {
			// We have intersection
			return t;
		}
	}
	return DBL_MAX;
}

void linearInterpolate(double xi, double yi, double zi, Triangle t, float p[3]) {
	double v1p[3], v2p[3], v3p[3];
	v1p[0] = t.v[0].position[0], v1p[1] = t.v[0].position[1], v1p[2] = t.v[0].position[2];
	v2p[0] = t.v[1].position[0], v2p[1] = t.v[1].position[1], v2p[2] = t.v[1].position[2];
	v3p[0] = t.v[2].position[0], v3p[1] = t.v[2].position[1], v3p[2] = t.v[2].position[2];
	double v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z;
	v1x = v1p[0] - xi, v1y = v1p[1] - yi, v1z = v1p[2] - zi;
	v2x = v2p[0] - xi, v2y = v2p[1] - yi, v2z = v2p[2] - zi;
	v3x = v3p[0] - xi, v3y = v3p[1] - yi, v3z = v3p[2] - zi;
	double v1Dist = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
	double v2Dist = sqrt(v2x * v2x + v2y * v2y + v2z * v2z);
	double v3Dist = sqrt(v3x * v3x + v3y * v3y + v3z * v3z);
	double totalDist = v1Dist + v2Dist + v3Dist;
	p[0] = (totalDist - v1Dist) / totalDist;
	p[1] = (totalDist - v2Dist) / totalDist;
	p[2] = (totalDist - v3Dist) / totalDist;
}

void draw_scene()
{
	unsigned int x,y;
	double increment = tan(fov / 2) * 2 / HEIGHT;
	double yInit = -tan(fov / 2);
	double xInit = WIDTH / HEIGHT * yInit;
	glPointSize(2.0);
	glBegin(GL_POINTS);
	// Get each ray <xVal, yVal, -1>
	for(x=0; x < WIDTH; x++)
	{
		double xVal = x * increment + xInit;
		for(y=0;y < HEIGHT;y++)
		{
			double yVal = y * increment + yInit;

			// Calculate normalized ray <ux, uy, uz>
			double magnitude = sqrt(xVal * xVal + yVal * yVal + 16);
			double xd = xVal / magnitude;
			double yd = yVal / magnitude;
			double zd = -4 / magnitude ;
			// printf("v = <%f,%f,%f>, |v| = %f\n", xd, yd, zd, sqrt(ux * ux + uy * uy + uz * uz));

			double minIntersectDistance = DBL_MAX;
			double minxi = 0, minyi = 0, minzi = 0;
			Triangle minTriangle;
			Sphere minSphere;
			bool minIsTri = false;

			// Check each ray intersection with triangles
			for (int i = 0; i < num_triangles; ++i) {
				Triangle triangle = triangles[i];
				float p[3], d[3];
				p[0] = 0, p[1] = 0, p[2] = 0;
				d[0] = xd, d[1] = yd, d[2] = zd;
				
				double intersect = rayIntersectsTriangle(p, d, triangle);
				if (intersect != DBL_MAX) {
					double xi = xd * intersect;
					double yi = yd * intersect;
					double zi = zd * intersect;
					double intersectDistance = xi * xi + yi * yi + zi * zi;
					if (intersectDistance < minIntersectDistance) {
						minIntersectDistance = intersectDistance;
						minxi = xi, minyi = yi, minzi = zi;
						minIsTri = true;
						minTriangle = triangle;
					}
				}
			}

			// Check each ray intersection with spheres
			for (int i = 0; i < num_spheres; i++) {
				Sphere sphere = spheres[i];
				float p[3], d[3];
				p[0] = 0, p[1] = 0, p[2] = 0;
				d[0] = xd, d[1] = yd, d[2] = zd;

				double intersect = rayIntersectsSphere(p, d, sphere);
				if (intersect != DBL_MAX) {
					double xi = xd * intersect;
					double yi = yd * intersect;
					double zi = zd * intersect;
					double intersectDistance = xi * xi + yi * yi + zi * zi;
					if (intersectDistance < minIntersectDistance) {
						minIntersectDistance = intersectDistance;
						minxi = xi, minyi = yi, minzi = zi;
						minIsTri = false;
						minSphere = sphere;
					}
				}
			}

			// Perform Phong shading on closest intersection
			if (minIntersectDistance < DBL_MAX) {
				float illumination[3];
				illumination[0] = ambient_light[0];
				illumination[1] = ambient_light[1];
				illumination[2] = ambient_light[2];
				// Send out shadow rays
				for (int i = 0; i < num_lights; ++i) {
					Light light = lights[i];
					float p[3], l[3];
					p[0] = minxi, p[1] = minyi, p[2] = minzi;
					l[0] = light.position[0] - minxi;
					l[1] = light.position[1] - minyi;
					l[2] = light.position[2] - minzi;

					// normalize shadow rays
					float magnitude = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
					l[0] = l[0] / magnitude;
					l[1] = l[1] / magnitude;
					l[2] = l[2] / magnitude;

					bool inShadow = false;
					for (int j = 0; j < num_triangles; ++j) {
						if (rayIntersectsTriangle(p, l, triangles[j]) < DBL_MAX) {
							inShadow = true;
							break;
						}
					}
					if (!inShadow) {
						for (int j = 0; j < num_spheres; ++j) {
							if (rayIntersectsSphere(p, l, spheres[j]) < DBL_MAX) {
								inShadow = true;
								break;
							}
						}
					}

					// Unblocked shadow ray
					if (!inShadow) {
						float n[3], v[3], r[3], temp[3], kd[3], ks[3];
						int sh = 0;
						if (minIsTri) {
							float tVec1[3], tVec2[3];
							vector(tVec1, minTriangle.v[0].position, p);
							vector(tVec2, minTriangle.v[1].position, p);
							crossProduct(n, tVec1, tVec2);
							magnitude = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
							n[0] = n[0] / magnitude;
							n[1] = n[1] / magnitude;
							n[2] = n[2] / magnitude;
							linearInterpolate(minxi, minyi, minzi, minTriangle, temp);
							kd[0] = minTriangle.v[0].color_diffuse[0] * temp[0] + minTriangle.v[1].color_diffuse[0] * temp[1] +
								minTriangle.v[2].color_diffuse[0] * temp[2];
							kd[1] = minTriangle.v[0].color_diffuse[1] * temp[0] + minTriangle.v[1].color_diffuse[1] * temp[1] +
								minTriangle.v[2].color_diffuse[1] * temp[2];
							kd[2] = minTriangle.v[0].color_diffuse[2] * temp[0] + minTriangle.v[1].color_diffuse[2] * temp[1] +
								minTriangle.v[2].color_diffuse[2] * temp[2];
							ks[0] = minTriangle.v[0].color_specular[0] * temp[0] + minTriangle.v[1].color_specular[0] * temp[1] +
								minTriangle.v[2].color_specular[0] * temp[2];
							ks[1] = minTriangle.v[0].color_specular[1] * temp[0] + minTriangle.v[1].color_specular[1] * temp[1] +
								minTriangle.v[2].color_specular[1] * temp[2];
							ks[2] = minTriangle.v[0].color_specular[2] * temp[0] + minTriangle.v[1].color_specular[2] * temp[1] +
								minTriangle.v[2].color_specular[2] * temp[2];
							sh = minTriangle.v[0].shininess * temp[0] + minTriangle.v[1].shininess * temp[1] + minTriangle.v[2].shininess * temp[2];
						}
						else {
							// Sphere normal, negate if inside sphere
							n[0] = (minxi - minSphere.position[0]) / minSphere.radius;
							n[1] = (minyi - minSphere.position[1]) / minSphere.radius;
							n[2] = (minzi - minSphere.position[2]) / minSphere.radius;
							kd[0] = minSphere.color_diffuse[0];
							kd[1] = minSphere.color_diffuse[1];
							kd[2] = minSphere.color_diffuse[2];
							ks[0] = minSphere.color_specular[0];
							ks[1] = minSphere.color_specular[1];
							ks[2] = minSphere.color_specular[2];
							sh = minSphere.shininess;
						}
						//printf("%f,%f,%f. %f\n", n[0], n[1], n[2], n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
						v[0] = -xd, v[1] = -yd, v[2] = -zd;
						scalarMultiply(temp, n, 2 * dotProduct(l, n));
						r[0] = temp[0] - l[0];
						r[1] = temp[1] - l[1];
						r[2] = temp[2] - l[2];
						// I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ sh)
						float lDotN = dotProduct(l, n);
						float rDotV = dotProduct(r, v);
						if (lDotN < 0)
							lDotN = 0;
						if (rDotV < 0)
							rDotV = 0;
						rDotV = pow(rDotV, sh);
						scalarMultiply(ks, ks, rDotV);
						scalarMultiply(kd, kd, lDotN);
						illumination[0] += light.color[0] * (kd[0] + ks[0]);
						illumination[1] += light.color[1] * (kd[1] + ks[1]);
						illumination[2] += light.color[2] * (kd[2] + ks[2]);
					}
				}
				if (illumination[0] > 1.0)
					illumination[0] = 1.0;
				if (illumination[1] > 1.0)
					illumination[1] = 1.0;
				if (illumination[2] > 1.0)
					illumination[2] = 1.0;
				//printf("%d,%d: (%f, %f, %f)\n", x, y, illumination[0] * 255, illumination[1] * 255, illumination[2] * 255);
				plot_pixel(WIDTH - x, HEIGHT - y, illumination[0] * 255, illumination[1] * 255, illumination[2] * 255);
			}
		}
	}
	glEnd();
	glFlush();
	printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

/* Write a jpg image from buffer*/
void save_jpg()
{
	if (filename == NULL)
		return;

	// Allocate a picture buffer // 
	cv::Mat3b bufferBGR = cv::Mat::zeros(HEIGHT, WIDTH, CV_8UC3); //rows, cols, 3-channel 8-bit.
	printf("File to save to: %s\n", filename);

	// unsigned char buffer[HEIGHT][WIDTH][3];
	for (int r = 0; r < HEIGHT; r++) {
		for (int c = 0; c < WIDTH; c++) {
			for (int chan = 0; chan < 3; chan++) {
				unsigned char red = buffer[r][c][0];
				unsigned char green = buffer[r][c][1];
				unsigned char blue = buffer[r][c][2];
				bufferBGR.at<cv::Vec3b>(r,c) = cv::Vec3b(blue, green, red);
			}
		}
	}
	if (cv::imwrite(filename, bufferBGR)) {
		printf("File saved Successfully\n");
	}
	else {
		printf("Error in Saving\n");
	}
}

void parse_check(char *expected,char *found)
{
  if(stricmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(stricmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(stricmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(stricmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{
	
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(1.0,1.0,1.0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      if(mode == MODE_JPEG)
	save_jpg();
    }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;
  
  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
