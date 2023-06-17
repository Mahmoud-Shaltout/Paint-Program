#if defined(UNICODE) && !defined(_UNICODE)
    #define _UNICODE
#elif defined(_UNICODE) && !defined(UNICODE)
    #define UNICODE
#endif


#include <vector>
#include <fstream>
#include <tchar.h>
#include <windows.h>
#include <math.h>
#include <iostream>
#include <vector>
#define MAXENTRIES 600
#include<climits>
#include <list>
#include <stack>

using namespace std;

#define File_save                    1

#define File_load                    2
#define File_choseColor              3
#define BK_color_white               4
#define BK_color_gray                5
#define BK_color_black               6
#define File_clear                   7
#define File_exit                    8

#define DrawLine_dda                 9
#define DrawLine_midline             10
#define DrawLine_paramline           11

#define DrawCircle_Direct            12
#define DrawCircle_Polar             13
#define DrawCircle_ItPolar           14
#define DrawCircle_midcircle         15
#define DrawCircle_midcircleModified 16

#define DrawEllipse_Direct           17
#define DrawEllipse_Polar            18
#define DrawEllipse_midellipse       19

#define Clipp_RecPoint               20
#define Clipp_RecLine                21
#define Clipp_RecPolygon             22
#define Clipp_SquPoint               23
#define Clipp_SquLine                24
#define Clipp_CirclePoint            25
#define Clipp_CircleLine             26

#define FillCircle_wline             27
#define FillCircle_wcircle           28
#define FillSquare                   29
#define FillRectangle                30
#define FillConvex                   31
#define FillnonConvex                32
#define FillRecursive                33
#define FillnonRecursive             34

#define square                       35
#define rectangle                    36

#define SplineCurve                  37


/* static variables */
static int x3 ;
static int y3 ;
POINT* pointArr ;
POINT dots [1000];
static int counter = 0 ;

vector<vector<int>> data;


/*  Declare Windows procedure  */
LRESULT CALLBACK WindowProcedure (HWND, UINT, WPARAM, LPARAM);

/*  Make the class name into a global variable  */
TCHAR szClassName[ ] = _T("CodeBlocksWindowsApp");

int WINAPI WinMain (HINSTANCE hThisInstance,
                     HINSTANCE hPrevInstance,
                     LPSTR lpszArgument,
                     int nCmdShow)
{
    HWND hwnd;               /* This is the handle for our window */
    MSG messages;            /* Here messages to the application are saved */
    WNDCLASSEX wincl;        /* Data structure for the windowclass */

    /* The Window structure */
    wincl.hInstance = hThisInstance;
    wincl.lpszClassName = szClassName;
    wincl.lpfnWndProc = WindowProcedure;      /* This function is called by windows */
    wincl.style = CS_DBLCLKS;                 /* Catch double-clicks */
    wincl.cbSize = sizeof (WNDCLASSEX);

    /* Use default icon and mouse-pointer */
    wincl.hIcon = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hIconSm = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hCursor = LoadCursor (NULL, IDC_CROSS);
    wincl.lpszMenuName = NULL;                 /* No menu */
    wincl.cbClsExtra = 0;                      /* No extra bytes after the window class */
    wincl.cbWndExtra = 0;                      /* structure or the window instance */
    /* Use Windows's default colour as the background of the window */
    wincl.hbrBackground = (HBRUSH) COLOR_BACKGROUND;

    /* Register the window class, and if it fails quit the program */
    if (!RegisterClassEx (&wincl))
        return 0;

    /* The class is registered, let's create the program*/
    hwnd = CreateWindowEx (
           0,                   /* Extended possibilites for variation */
           szClassName,         /* Classname */
           _T("Code::Blocks Template Windows App"),       /* Title Text */
           WS_OVERLAPPEDWINDOW, /* default window */
           CW_USEDEFAULT,       /* Windows decides the position */
           CW_USEDEFAULT,       /* where the window ends up on the screen */
           700,                 /* The programs width 544 default value*/
           700,                 /* and height in pixels 375 default value*/
           HWND_DESKTOP,        /* The window is a child-window to desktop */
           NULL,                /* No menu */
           hThisInstance,       /* Program Instance handler */
           NULL                 /* No Window Creation data */
           );

    /* Make the window visible on the screen */
    ShowWindow (hwnd, nCmdShow);

    /* Run the message loop. It will run until GetMessage() returns 0 */
    while (GetMessage (&messages, NULL, 0, 0))
    {
        /* Translate virtual-key messages into character messages */
        TranslateMessage(&messages);
        /* Send message to WindowProcedure */
        DispatchMessage(&messages);
    }

    /* The program return-value is 0 - The value that PostQuitMessage() gave */
    return messages.wParam;
}

void swap(int& x1, int& y1f, int& x2, int& y2)
{
	int tmp = x1;
	x1 = x2;
	x2 = tmp;
	tmp = y1f;
	y1f = y2;
	y2 = tmp;
}
int Round(double x)
{
	return (int)(x + 0.5);
}
int calcR(int x1, int y1, int x2, int y2){
    int dx = x2 - x1;
    int dy = y2 - y1;

    return sqrt((dx*dx) + (dy*dy));
}
COLORREF color(int colorChoice){
    if(colorChoice == 1)
        return RGB(0,0,0);
    else if(colorChoice == 2)
        return RGB(255,0,0);
    else if(colorChoice == 3)
        return RGB(0,255,0);
    else if(colorChoice == 4)
        return RGB(0,0,255);
    else if(colorChoice == 5)
        return RGB(255,255,0);
    else if(colorChoice == 6)
        return RGB(255,0,255);
    else if(colorChoice == 7)
        return RGB(0,255,255);
    else if(colorChoice == 8)
        return RGB(255,255,255);
}


struct point
{
    int x;
    int y;
};

struct shape
{
    int type;
    vector<point> data;
    COLORREF color;
};

point current_point;
shape current_shape;
vector<shape> shapes;
void get_RGB_values(COLORREF color, int& R, int& G, int& B)
{
    R = color & 0xff;
    G = (color >> 8) & 0xff;
    B = (color >> 16) & 0xff;
}
void print_shapes(){ // for debugging only
    int r, g, b;
    for (int i = 0; i < shapes.size(); ++i)
    {
        cout << "TYPE = " << shapes[i].type << endl;
        cout << "POINTS = ";
        for (int j = 0; j < shapes[i].data.size(); ++j)
            cout << "(" << shapes[i].data[j].x << ", " << shapes[i].data[j].y << ")" << " ";

        cout << endl;
        get_RGB_values(shapes[i].color, r, g, b);
        cout << "COLOR = ";
        cout << r << " ";
        cout << g << " ";
        cout << b << endl;
        cout << "============" << endl;
    }
}

struct Vertex
{
	double x, y;
	Vertex(int x1 = 0, int y1 = 0)
	{
		x = x1;
		y = y1;
	}
};
typedef vector<Vertex> VertexList;
typedef bool (*IsInFunc)(Vertex& v, int edge);
typedef Vertex(*IntersectFunc)(Vertex& v1, Vertex& v2, int edge);



/*drawing line*/
//Line DDA
void DrawLineDDA(HDC hdc, int x1, int y1f, int x2, int y2, COLORREF c)
{
	int dx = x2 - x1;
	int dy = y2 - y1f;
	if (abs(dy) <= abs(dx))
	{
		if (x1 > x2)swap(x1, y1f, x2, y2);
        SetPixel(hdc, x1, y1f, c);
		int x = x1;
		while (x < x2)
		{
			x++;
			double y = y1f + (double)(x - x1)*dy / dx;
			SetPixel(hdc, x, Round(y), c);
		}
	}
	else {
		if (y1f > y2)swap(x1, y1f, x2, y2);
		SetPixel(hdc, x1, y1f, c);
		int y = y1f;
		while (y < y2)
		{
			y++;
			double x = x1 + (double)(y - y1f)*dx / dy;
			SetPixel(hdc, Round(x), y, c);
		}
	}

}
//Line Midpoint
void DrawLineMidpoint(HDC hdc, int x1, int y1f, int x2, int y2, COLORREF c)
{
    int x = x1;
    int y = y1f;
    int dx = x2-x1;
    int dy = y2-y1f;
    int x_dir = 1;
    int y_dir = 1;
    if(dx<0){
        dx = -dx;
        x_dir = -1;
    }
    if(dy<0){
        dy = -dy;
        y_dir = -1;
    }
    if(dy>dx){
        swap(x1,x2,y1f,y2);
    }
    if(dx>dy){
        int p =2*dy-dx;
        int two_dy = 2*dy;
        int two_minus= 2*dy - 2*dx;
        for(int i=1; i<=dx; i++){
            x+=x_dir;
            if(p<0)
                p += two_dy;
            else{
                p += two_minus;
                y += y_dir;
            }
            SetPixel(hdc,x,y,c);
        }
    }else{
        int p =2*dx-dy;
        int two_dy = 2*dx;
        int two_minus= 2*dx - 2*dy;
        for(int i=1; i<=dy; i++){
            y+=y_dir;
            if(p<0)
                p += two_dy;
            else{
                p += two_minus;
                x += x_dir;
            }
            SetPixel(hdc,x,y,c);
        }
    }
}
//Line Parametric
void DrawLineParametric(HDC hdc, int x1, int y1f, int x2, int y2, COLORREF c)
{
    double dx = x2 - x1;
	double dy = y2 - y1f;
	int x;
	int y;
	for(double t=0;t<1;t = t + 0.001){
        x = x1 + t * dx;
		y = y1f + t * dy;
		SetPixel(hdc,Round(x) , Round(y) ,c );
	}
}

/*drawing circle*/
//8Points
void Draw8Points(HDC hdc, int xc, int yc, int x, int y, COLORREF c)
{
	SetPixel(hdc, xc + x, yc + y, c);
	SetPixel(hdc, xc + x, yc - y, c);
	SetPixel(hdc, xc - x, yc - y, c);
	SetPixel(hdc, xc - x, yc + y, c);
	SetPixel(hdc, xc + y, yc + x, c);
	SetPixel(hdc, xc + y, yc - x, c);
	SetPixel(hdc, xc - y, yc - x, c);
	SetPixel(hdc, xc - y, yc + x, c);


}
void copyPoints ()
{
    pointArr = new POINT [counter];
    for (int i = 0 ; i < counter ; i++)
    {
        pointArr[i].x = dots[i].x ;
        pointArr[i].y = dots[i].y ;
    }
}
//CircleDirect
void DrawCircleDirect(HDC hdc, int xc, int yc, int R, COLORREF c,int quarter)
{
	int x = 0,y = R;
	while (x < y)
	{
		x++;
		y = sqrt(R * R - x * x);
		if(quarter==0){
            Draw8Points(hdc, xc, yc, 0, R, c);
            Draw8Points(hdc, xc, yc, x, Round(y), c);
		}else if(quarter==1){
            SetPixel(hdc, xc + x, yc - y, c);
            SetPixel(hdc, xc + y, yc - x, c);
		}else if(quarter==2){
            SetPixel(hdc, xc - x, yc - y, c);
            SetPixel(hdc, xc - y, yc - x, c);
		}else if(quarter==3){
            SetPixel(hdc, xc - x, yc + y, c);
            SetPixel(hdc, xc - y, yc + x, c);
		}else if(quarter==4){
		    SetPixel(hdc, xc + x, yc + y, c);
            SetPixel(hdc, xc + y, yc + x, c);
		}
	}
}
//CircleMidPoint
void DrawCircleMidPoint(HDC hdc, int xc, int yc, int R, COLORREF c,int quarter)
{
    double x = 0 ;
    double y = R ;
    Draw8Points(hdc,xc, yc, x, y,c);
    while(x<y)
    {
        double d = (x+1)*(x+1) + (y-0.5) * (y-0.5) - R*R;
        if (d < 0 )
        {
            x++;
        }
        else
        {
            x++;
            y--;
        }
        Draw8Points(hdc,xc, yc, x, y,c);
    }

}
//CircleMidPoint_modified
void DrawCircleMidPoint_modified(HDC hdc, int xc, int yc, int R, COLORREF c,int quarter)
{
	double x = 0;
	double y = R;
	double d = 1 - R;
	Draw8Points(hdc,xc, yc, x, y,c);
	while (x < y)
	{
		if (d < 0)
		{
			d = d + 2 * x + 3;
			x++;
		}
		else
		{
			d = d + 2 * (x - y) + 5;
			x++;
			y--;
		}
		Draw8Points(hdc,xc, yc, x, y,c);
	}
}
//CirclePolar
void DrawCirclePolar(HDC hdc, int xc, int yc, int R, COLORREF c,int quarter)
{
    int x=R;
    int y=0;
    double theta=0;
    double dtheta=1.0/R;
    Draw8Points(hdc,xc,yc,x,y,c);
    while(x>y)
    {
        theta+=dtheta;
        x=round(R*cos(theta));
        y=round(R*sin(theta));
        Draw8Points(hdc,xc,yc,x,y,c);
    }
}
//CircleIterativePolar
void DrawCircleIterativePolar(HDC hdc,int xc,int yc, int R,COLORREF c,int quarter)
{
    double x=R,y=0;
    double dtheta=1.0/R;
    double cdtheta=cos(dtheta),sdtheta=sin(dtheta);
    Draw8Points(hdc,xc,yc,R,0,c);
    while(x>y)
    {
        double x1=x*cdtheta-y*sdtheta;
        y=x*sdtheta+y*cdtheta;
        x=x1;
        Draw8Points(hdc,xc,yc,round(x),round(y),c);
    }
}

/*drawing ellipse*/
//4Points
void Draw4Points(HDC hdc, int xc, int yc, int x, int y, COLORREF c)
{
	SetPixel(hdc, xc + x, yc + y, c);
	SetPixel(hdc, xc + x, yc - y, c);
	SetPixel(hdc, xc - x, yc + y, c);
	SetPixel(hdc, xc - x, yc - y, c);
}
//Ellipse Direct
void DrawEllipseDirect(HDC hdc, int xc, int yc, int A, int B, COLORREF c) {
	int x = 0;
	double y = B;
	Draw4Points(hdc, xc, yc, 0, B, c);
	while (x * B * B < y * A * A) {
		x++;
		y = B * sqrt(1.0 - (double)x * x / (A * A));
		Draw4Points(hdc, xc, yc, x, Round(y), c);
	}
	int y1 = 0;
	double x1 = A;
	Draw4Points(hdc, xc, yc, A, 0, c);

	while (y1 * A * A < x1 * B * B) {
		y1++;
		x1 = A * sqrt(1.0 - (double)y1 * y1 / (B * B));
		Draw4Points(hdc, xc, yc, Round(x1), y1, c);
	}
}
//Ellipse Polar
void DrawEllipsePolar(HDC hdc, int xc, int yc, int A, int B, COLORREF c) {
	double x = 0;
	double y = B;
	double theta = 0;
	double range = A;
	Draw4Points(hdc, xc, yc, 0, B, c);
	while (theta <= range) {
		theta += 1.0 / max(A, B);
		x = A * cos(theta);
		y = B * sin(theta);
		Draw4Points(hdc, xc, yc, Round(x), Round(y), c);
	}
}
// Ellipse mid_point
void midellipse(HDC hdc,float xc,float yc,float Rx,float Ry,COLORREF c){
    float x,y;

    float P1 = pow(Ry, 2) - pow(Rx, 2)* Ry + (1/4)*pow(Rx,2);
        x = 0;
        y = Ry;
        float dx=2*pow(Ry,2)*x;
        float dy=2*pow(Rx,2)*y;
        while (dx < dy)
        {
            if( P1 < 0)
            {
                x++;
                dx=dx+(2 * Ry * Ry);
                P1 += 2 * pow(Ry,2)* x + pow(Ry,2);
            }
            else
            {
                x++;
                y--;
                dx = dx + (2 * Ry * Ry);
                dy = dy - (2 * Rx * Rx);
                P1 += 2 * pow(Ry,2)* x - 2 * pow(Rx, 2)* y + pow(Ry, 2);
            }
            Draw4Points(hdc,xc,yc,x,y,c);

        }

        float P2 = pow(Ry, 2) * pow(x + 0.5,2) + pow(Rx, 2) * pow(y-0.5, 2) - pow(Rx,2) * pow(Ry, 2);
        while (y >= 0)
        {
            if (P2 > 0)
            {
                y--;
                dy = dy - (2 * Rx * Rx);
                P2 += pow(Rx,2) - 2*pow(Rx, 2)*y;
            }
            else
            {
                x++;
                y--;
                 dx = dx + (2 * Ry * Ry);
                 dy = dy - (2 * Rx * Rx);
                P2 += 2*pow(Ry, 2)*x - 2*pow(Rx, 2)*y + pow(Rx,2);
            }
            Draw4Points(hdc,xc,yc,x,y,c);
        }

}

/*clipping*/

//point rectangle clipping
union OutCode
{
	unsigned All : 4;
	struct { unsigned left : 1,top : 1,right : 1,bottom : 1; };
};
OutCode GetOutCode(double x, double y, int xleft, int ytop, int xright, int ybottom)
{
	OutCode out;
	out.All = 0;
	if (x < xleft)out.left = 1; else if (x > xright)out.right = 1;
	if (y < ytop)out.top = 1; else if (y > ybottom)out.bottom = 1;
	return out;
}
void VIntersect(double xs, double ys, double xe, double ye, int x, double* xi, double* yi)
{
	*xi = x;
	*yi = ys + (x - xs) * (ye - ys) / (xe - xs);
}
void HIntersect(double xs, double ys, double xe, double ye, int y, double* xi, double* yi)
{
	*yi = y;
	*xi = xs + (y - ys) * (xe - xs) / (ye - ys);
}
void CohenSuth(HDC hdc, int xs, int ys, int xe, int ye, int xleft, int ytop, int xright, int ybottom)
{
	double x1 = xs, y1 = ys, x2 = xe, y2 = ye;
	OutCode out1 = GetOutCode(x1, y1, xleft, ytop, xright, ybottom);
	OutCode out2 = GetOutCode(x2, y2, xleft, ytop, xright, ybottom);
	while ((out1.All || out2.All) && !(out1.All & out2.All))
	{
		double xi, yi;
		if (out1.All)
		{
			if (out1.left)VIntersect(x1, y1, x2, y2, xleft, &xi, &yi);
			else if (out1.top)HIntersect(x1, y1, x2, y2, ytop, &xi, &yi);
			else if (out1.right)VIntersect(x1, y1, x2, y2, xright, &xi, &yi);
			else HIntersect(x1, y1, x2, y2, ybottom, &xi, &yi);
			x1 = xi;
			y1 = yi;
			out1 = GetOutCode(x1, y1, xleft, ytop, xright, ybottom);
		}
		else
		{
			if (out2.left)VIntersect(x1, y1, x2, y2, xleft, &xi, &yi);
			else if (out2.top)HIntersect(x1, y1, x2, y2, ytop, &xi, &yi);
			else if (out2.right)VIntersect(x1, y1, x2, y2, xright, &xi, &yi);
			else HIntersect(x1, y1, x2, y2, ybottom, &xi, &yi);
			x2 = xi;
			y2 = yi;
			out2 = GetOutCode(x2, y2, xleft, ytop, xright, ybottom);
		}
	}
	if (!out1.All && !out2.All)
	{
		MoveToEx(hdc, Round(x1), Round(y1), NULL);
		LineTo(hdc, Round(x2), Round(y2));
	}
}
void PointClipping(HDC hdc, int x1, int y1, int X_left, int Y_top, int X_right, int Y_bottom) {
	if (x1 >= X_left && x1 <= X_right && y1 >= Y_top && y1 <= Y_bottom) {
		SetPixel(hdc, x1, y1, RGB(0, 0, 0));
	}
}



//point circle clipping
void PointClippingInCircle(HDC hdc, int x1, int y1f, int xc, int yc, int r1)
{
	int x2 = (x1 - xc) * (x1 - xc);
	int y2 = (y1f - yc) * (y1f - yc);
	int R2 = r1 * r1;
	double d = x2 + y2 - R2;
	if (d <= 0) {
		SetPixel(hdc, x1, y1f, RGB(0, 0, 0));
	}
}
//line circle clipping
void LineClippingInCircle(HDC hdc, int x1, int y1f, int x2, int y2, int xc, int yc, int r1)
{
	COLORREF c = RGB(0, 0, 0);
	int dx = x2 - x1;
	int dy = y2 - y1f;
	if (abs(dy) <= abs(dx))
	{
		if (x1 > x2) swap(x1, y1f, x2, y2);
		PointClippingInCircle(hdc, x1, y1f, xc, yc, r1);
		int x = x1;
		while (x < x2)
		{
			x++;
			double y = y1f + (double)(x - x1) * dy / dx;
			PointClippingInCircle(hdc, x, y, xc, yc, r1);
		}
	}
	else {
		if (y1f > y2)swap(x1, y1f, x2, y2);
		PointClippingInCircle(hdc, x1, y1f, xc, yc, r1);
		int y = y1f;
		while (y < y2)
		{
			y++;
			double x = x1 + (double)(y - y1f) * dx / dy;
			PointClippingInCircle(hdc, x, y, xc, yc, r1);
		}
	}
}

/*fill circle */
void fillingWithLines(HDC hdc, int xc, int yc, int R, COLORREF c,int quarter)
{
	int x = 0,y = R;
	while (x < y)
	{
		x++;
		y = sqrt(R * R - x * x);
		if(quarter==1){
            DrawLineDDA(hdc, xc + y, yc - Round(x), xc, yc, c);
            DrawLineDDA(hdc, xc + x, yc - Round(y), xc, yc, c);
		}else if(quarter==2){
            DrawLineDDA(hdc, xc - y, yc - Round(x), xc, yc, c);
            DrawLineDDA(hdc, xc - x, yc - Round(y), xc, yc, c);
		}else if(quarter==3){
            DrawLineDDA(hdc, xc - y, yc + Round(x), xc, yc, c);
            DrawLineDDA(hdc, xc - x, yc + Round(y), xc, yc, c);
		}else if(quarter==4){
            DrawLineDDA(hdc, xc + y, yc + Round(x), xc, yc, c);
            DrawLineDDA(hdc, xc + x, yc + Round(y), xc, yc, c);
		}else
            cout<< "Invalid Input!!";

	}
}
void fillingWithCircles(HDC hdc, int xc, int yc, int R, COLORREF c,int quarter)
{
	while (R>0)
	{
		DrawCircleDirect(hdc, xc, yc, R--, c,quarter);
	}
}

VertexList ClipWithEdge(VertexList p, int edge, IsInFunc In, IntersectFunc Intersect)
{
	VertexList OutList;
	Vertex v1 = p[p.size() - 1];
	bool v1_in = In(v1, edge);
	for (int i = 0; i < (int)p.size(); i++)
	{
		Vertex v2 = p[i];
		bool v2_in = In(v2, edge);
		if (!v1_in && v2_in)
		{
			OutList.push_back(Intersect(v1, v2, edge));
			OutList.push_back(v2);
		}
		else if (v1_in && v2_in) OutList.push_back(v2);
		else if (v1_in) OutList.push_back(Intersect(v1, v2, edge));
		v1 = v2;
		v1_in = v2_in;
	}
	return OutList;
}
bool InLeft(Vertex& v, int edge)
{
	return v.x >= edge;
}
bool InRight(Vertex& v, int edge)
{
	return v.x <= edge;
}
bool InTop(Vertex& v, int edge)
{
	return v.y >= edge;
}
bool InBottom(Vertex& v, int edge)
{
	return v.y <= edge;
}
Vertex VIntersect(Vertex& v1, Vertex& v2, int xedge)
{
	Vertex res;
	res.x = xedge;
	res.y = v1.y + (xedge - v1.x) * (v2.y - v1.y) / (v2.x - v1.x);
	return res;
}
Vertex HIntersect(Vertex& v1, Vertex& v2, int yedge)
{
	Vertex res;
	res.y = yedge;
	res.x = v1.x + (yedge - v1.y) * (v2.x - v1.x) / (v2.y - v1.y);
	return res;
}
void PolygonClip(HDC hdc, Vertex p[], int n, int xleft, int ytop, int xright, int ybottom)
{
	VertexList vlist;
	for (int i = 0; i < n; i++)vlist.push_back(Vertex(p[i].x, p[i].y));
	vlist = ClipWithEdge(vlist, xleft, InLeft, VIntersect);
	vlist = ClipWithEdge(vlist, ytop, InTop, HIntersect);
	vlist = ClipWithEdge(vlist, xright, InRight, VIntersect);
	vlist = ClipWithEdge(vlist, ybottom, InBottom, HIntersect);
	Vertex v1 = vlist[vlist.size() - 1];
	for (int i = 0; i < (int)vlist.size(); i++)
	{
		Vertex v2 = vlist[i];
		MoveToEx(hdc, Round(v1.x), Round(v1.y), NULL);
		LineTo(hdc, Round(v2.x), Round(v2.y));
		v1 = v2;
	}
}


struct Entry
{
    int xmin,xmax;
};
void InitEntries(Entry table[])
{
    for(int i = 0; i < MAXENTRIES ; i++ )
    {
    table[i].xmin =  INT_MAX;
    table[i].xmax = -INT_MAX;
    }
}

void ScanEdge(POINT v1,POINT v2,Entry table[])
{
    if(v1.y==v2.y)return;
    if(v1.y>v2.y)swap(v1,v2);
    double minv=(double)(v2.x-v1.x)/(v2.y-v1.y);
    double x=v1.x;
    int y=v1.y;
    while(y<v2.y)
    {
        if(x<table[y].xmin)table[y].xmin=(int)ceil(x);
        if(x>table[y].xmax)table[y].xmax=(int)floor(x);
        y++;
        x+=minv;
    }
}
void DrawSanLines(HDC hdc,Entry table[],COLORREF color)
{
for(int y=0;y<MAXENTRIES;y++)
    if(table[y].xmin<table[y].xmax)
        for(int x=table[y].xmin;x<=table[y].xmax;x++)
            SetPixel(hdc,x,y,color);
}

// convex filling
void ConvexFill(HDC hdc,POINT p[],int n,COLORREF color)
{
    Entry *table=new Entry[MAXENTRIES];
    InitEntries(table);
    POINT v1=p[n-1];
    for(int i=0;i<n;i++)
        {
            POINT v2=p[i];
            ScanEdge(v1,v2,table);
            v1=p[i];
        }
DrawSanLines(hdc,table,color);
delete table;
}

// general Polygon Filling :
// some utilities :
struct EdgeRec
{
    double x;
    double minv;
    int ymax;
    bool operator<(EdgeRec r)
    {
        return x<r.x;
    }
};
typedef list<EdgeRec> EdgeList;

EdgeRec InitEdgeRec(POINT& v1,POINT& v2)
{
    if(v1.y>v2.y)swap(v1,v2);
    EdgeRec rec;
    rec.x=v1.x;
    rec.ymax=v2.y;
    rec.minv=(double)(v2.x-v1.x)/(v2.y-v1.y);
return rec;
}

void InitEdgeTable(POINT *polygon,int n,EdgeList table[])
{
POINT v1=polygon[n-1];
    for(int i=0;i<n;i++)
    {
        POINT v2=polygon[i];
        if(v1.y==v2.y){v1=v2;continue;}
        EdgeRec rec=InitEdgeRec(v1, v2);
        table[v1.y].push_back(rec);
        v1=polygon[i];
    }
}
// the main algo :
void GeneralPolygonFill(HDC hdc,POINT *polygon,int n,COLORREF c)
{
    EdgeList *table=new EdgeList [MAXENTRIES];
    InitEdgeTable(polygon,n,table);
    int y=0;

    while(y<MAXENTRIES && table[y].size()==0)y++;

    if(y==MAXENTRIES)return;
    EdgeList ActiveList=table[y];
    while (ActiveList.size()>0)
    {
        ActiveList.sort();
        for(EdgeList::iterator it=ActiveList.begin();it!=ActiveList.end();it++)
        {
            int x1=(int)ceil(it->x);
            it++;
            int x2=(int)floor(it->x);
            for(int x=x1;x<=x2;x++)SetPixel(hdc,x,y,c);
        }
        y++;
        EdgeList::iterator it=ActiveList.begin();
        while(it!=ActiveList.end())
        if(y==it->ymax) it=ActiveList.erase(it); else it++;
        for(EdgeList::iterator it=ActiveList.begin();it!=ActiveList.end();it++)
        it->x+=it->minv;
        ActiveList.insert(ActiveList.end(),table[y].begin(),table[y].end());
    }
    delete[] table;
}



// recursive flood fill :
void FloodFill(HDC hdc,int x,int y,COLORREF Cb,COLORREF Cf)
{
    COLORREF C=GetPixel(hdc,x,y);

    if(C==Cb || C==Cf)return;

    SetPixel(hdc,x,y,Cf);
    FloodFill(hdc,x+1,y,Cb,Cf);
    FloodFill(hdc,x-1,y,Cb,Cf);
    FloodFill(hdc,x,y+1,Cb,Cf);
    FloodFill(hdc,x,y-1,Cb,Cf);
}

// non-recursive flood fill :
void NRFloodFill(HDC hdc,int x,int y,COLORREF Cb,COLORREF Cf)
{
    stack<Vertex> S;
    S.push(Vertex(x,y));
    while(!S.empty())
    {
        Vertex v=S.top();
        S.pop();

        COLORREF c=GetPixel(hdc,v.x,v.y);

        if(c==Cb || c==Cf)continue;

        SetPixel(hdc,v.x,v.y,Cf);
        S.push(Vertex(v.x+1,v.y));
        S.push(Vertex(v.x-1,v.y));
        S.push(Vertex(v.x,v.y+1));
        S.push(Vertex(v.x,v.y-1));
    }
}


void AddMenus(HWND);
void BKcolor(HWND hwnd,COLORREF BKc);
HMENU hmenu;

int X_left, Y_top, X_right, Y_bottom;


//COLORREF c = RGB(255, 0, 0);
COLORREF BKc;
int colorChoice,quarter=0;
int choice=0;
int t2;

// Saving and Loading To-From File
bool HDCToFile(const char* FilePath, HDC Context, RECT Area, uint16_t BitsPerPixel = 24)
{
    uint32_t Width = Area.right - Area.left;
    uint32_t Height = Area.bottom - Area.top;
    BITMAPINFO Info;
    BITMAPFILEHEADER Header;
    memset(&Info, 0, sizeof(Info));
    memset(&Header, 0, sizeof(Header));
    Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    Info.bmiHeader.biWidth = Width;
    Info.bmiHeader.biHeight = Height;
    Info.bmiHeader.biPlanes = 1;
    Info.bmiHeader.biBitCount = BitsPerPixel;
    Info.bmiHeader.biCompression = BI_RGB;
    Info.bmiHeader.biSizeImage = Width * Height * (BitsPerPixel > 24 ? 4 : 3);
    Header.bfType = 0x4D42;
    Header.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
    char* Pixels = NULL;
    HDC MemDC = CreateCompatibleDC(Context);
    HBITMAP Section = CreateDIBSection(Context, &Info, DIB_RGB_COLORS, (void**)&Pixels, 0, 0);
    DeleteObject(SelectObject(MemDC, Section));
    BitBlt(MemDC, 0, 0, Width, Height, Context, Area.left, Area.top, SRCCOPY);
    DeleteDC(MemDC);
    std::fstream hFile(FilePath, std::ios::out | std::ios::binary);
    if (hFile.is_open())
    {
        hFile.write((char*)&Header, sizeof(Header));
        hFile.write((char*)&Info.bmiHeader, sizeof(Info.bmiHeader));
        hFile.write(Pixels, (((BitsPerPixel * Width + 31) & ~31) / 8) * Height);
        hFile.close();
        DeleteObject(Section);
        return true;
    }
    DeleteObject(Section);
    return false;
}
void load(HWND hWnd, HDC &hdc)
{
    char fileName[12] = "picture.bmp";
    if (fileName == "")
        return ;
    HBITMAP hBitmap;
    hBitmap = (HBITMAP)::LoadImage(NULL, fileName, IMAGE_BITMAP, 0, 0, LR_LOADFROMFILE);
    HDC hLocalDC;
    hLocalDC = CreateCompatibleDC(hdc);
    BITMAP qBitmap;
    int iReturn = GetObject(reinterpret_cast<HGDIOBJ>(hBitmap), sizeof(BITMAP),reinterpret_cast<LPVOID>(&qBitmap));
    HBITMAP hOldBmp = (HBITMAP)SelectObject(hLocalDC, hBitmap);
    BOOL qRetBlit = BitBlt(hdc, 0, 0, qBitmap.bmWidth, qBitmap.bmHeight,hLocalDC, 0, 0, SRCCOPY);
    SelectObject (hLocalDC, hOldBmp);
    DeleteDC(hLocalDC);
    DeleteObject(hBitmap);
}
void save(HWND &hWnd)
{
    HDC hdc = GetDC(hWnd);
    char fileName[12] = "picture.bmp";
    if (fileName == "")
        return ;
    int windowWidth ;
    int windowHeight;
    RECT rect;
    if(GetWindowRect(hWnd, &rect))
    {
        windowWidth = rect.right - rect.left;
        windowHeight = rect.bottom - rect.top;
    }
    RECT rect1 = {0, 0, windowWidth, windowHeight};
    HDCToFile(fileName,hdc,rect1);
}





static int B ;
static int A ;
static int t = 1, xc1, xc2, yc1, yc2, x, r1, r2;
static int x1, y1f, x2, y2;
static int xc,yc,y,R,yy1;


/* bezier curve */
struct Vector {
	double v[2];
	Vector(double x = 0, double y = 0)
	{
		v[0] = x; v[1] = y;
	}
	double& operator[](int i) {
		return v[i];
	}
};
void DrawBezierCurve(HDC hdc, Vector& p1, Vector& p2, Vector& p3, Vector& p4, COLORREF c)
{
	double a0 = p1[0], a1 = 3 * (p2[0] - p1[0] ),
		a2 = 3 * p1[0] - 6 * p2[0] + 3 * p3[0],
		a3 =   -p1[0] + 3 * p2[0] - 3 * p3[0] + p4[0];
	double b0 = p1[1], b1 = 3 * (p2[1] - p1[1]),
		b2 = 3 * p1[1] - 6 * p2[1] + 3 * p3[1],
		b3 = -p1[1] + 3 * p2[1] - 3 * p3[1] + p4[1];
	for (double t = 0; t <= 1; t += 0.001)
	{
		double t2 = t * t, t3 = t2 * t;
		double x = a0 + a1 * t + a2 * t2 + a3 * t3;
		double y = b0 + b1 * t + b2 * t2 + b3 * t3;
		SetPixel(hdc, Round(x), Round(y), c);
	}
}

/*  Hermite curve */
void DrawHermiteCurve(HDC hdc,Vector& p1, Vector& T1, Vector& p2, Vector& T2, COLORREF c)
{
	double a0 = p1[0], a1 = T1[0],
		a2 = -3 * p1[0] - 2 * T1[0] + 3 * p2[0] - T2[0],
		a3 = 2 * p1[0] + T1[0] - 2 * p2[0] + T2[0];
	double b0 = p1[1], b1 = T1[1],
		b2 = -3 * p1[1] - 2 * T1[1] + 3 * p2[1] - T2[1],
		b3 = 2 * p1[1] + T1[1] - 2 * p2[1] + T2[1];
	for (double t = 0; t <= 1; t += 0.001)
	{
		double t2 = t*t, t3 = t2*t;
		double x = a0 + a1*t + a2*t2 + a3*t3;
		double y = b0 + b1*t + b2*t2 + b3*t3;
		SetPixel(hdc, Round(x), Round(y), c);
	}
}

/*draw spline curve*/
void DrawCardinalSpline(HDC hdc,Vector P[],int n,double c,COLORREF color) {
    double c1=1-c;
    Vector T0(c1*(P[2].v[0]-P[0].v[0]),c1*(P[2].v[1]-P[0].v[1]));
    for(int i=2;i<n-1;i++){
        SetPixel(hdc,P[i].v[0],P[i].v[1],color);
        Vector T1(c1*(P[i+1].v[0]-P[i-1].v[0]),c1*(P[i+1].v[1]-P[i-1].v[1]));
        DrawHermiteCurve(hdc,P[i-1],T0,P[i],T1,color);
        T0=T1;
    }
}

/*  This function is called by the Windows function DispatchMessage()  */
LRESULT CALLBACK WindowProcedure (HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    static HDC hdc;
    static int B ;
    static int A ;
    static int t = 1, xc1, xc2, yc1, yc2, x, r1, r2;
    static int x1, y1f, x2, y2;
    static int xc,yc,y,R,yy1;
	static Vector p[4];
	static Vertex p1[20];
	static Vector splinepoint[6];

    switch (message)                  /* handle the messages */
    {
        case WM_COMMAND:
            switch(LOWORD(wParam))
            {
            /*file menu*/
            case File_save:
                save(hwnd);
                break;
            case File_load:
                load(hwnd,hdc);
                break;
            case File_choseColor:
                cout<< "Enter color number: \n1-Black \n2-Red \n3-Green \n4-Blue \n5-Yellow \n6-Purple \n7-Cyan \n8-White \n";
                cin>> colorChoice;
                break;
            case BK_color_white:
                cout<< "Background color changed to be white!!"<<endl;
                BKc = RGB(255,255,255);
                BKcolor(hwnd,BKc);
                break;
            case BK_color_gray:
                cout<< "Background color changed to be gray!!"<<endl;
                BKc = RGB(200,200,200);
                BKcolor(hwnd,BKc);
                break;
            case BK_color_black:
                cout<< "Background color changed to be black!!"<<endl;
                BKc = RGB(0,0,0);
                BKcolor(hwnd,BKc);
                break;
            case File_clear:
                InvalidateRect(hwnd, NULL, TRUE);
                break;
            case File_exit:
                char decide;
                cout << "Are you want to exit? (y/n)" << endl;
                cin >> decide;
                if(decide == 'y')
                    SendMessage(hwnd, WM_CLOSE, NULL, NULL);
                break;

            /*draw menu*/
            //draw line
            case DrawLine_dda:
                choice = DrawLine_dda;
                break;
            case DrawLine_midline:
                choice = DrawLine_midline;
                break;
            case DrawLine_paramline:
                choice = DrawLine_paramline;
                break;

            //draw circle
            case DrawCircle_Direct:
                choice = DrawCircle_Direct;
                break;
            case DrawCircle_Polar:
                choice = DrawCircle_Polar;
                break;
            case DrawCircle_ItPolar:
                choice = DrawCircle_ItPolar;
                break;
            case DrawCircle_midcircle:
                choice = DrawCircle_midcircle;
                break;
            case DrawCircle_midcircleModified:
                choice = DrawCircle_midcircleModified;
                break;

            //draw ellipse
            case DrawEllipse_Direct:
                choice = DrawEllipse_Direct;
                break;
            case DrawEllipse_Polar:
                choice = DrawEllipse_Polar;
                break;
            case DrawEllipse_midellipse:
                choice = DrawEllipse_midellipse;
                break;

            /*clipping*/
            //rectangle clipping
            case Clipp_RecPoint:
                choice = Clipp_RecPoint;
                break;
            case Clipp_RecLine:
                choice = Clipp_RecLine;
                break;
            case Clipp_RecPolygon:
                choice = Clipp_RecPolygon;
                break;
            //square clipping
            case Clipp_SquPoint:
                choice = Clipp_SquPoint;
                break;
            case Clipp_SquLine:
                choice = Clipp_SquLine;
                break;
            //circle clipping
            case Clipp_CirclePoint:
                choice = Clipp_CirclePoint;
                break;
            case Clipp_CircleLine:
                choice = Clipp_CircleLine;
                break;
            case SplineCurve:
                choice = SplineCurve;
                break;


            /*filling*/
            //circle filling
            case FillCircle_wline:
                cout<< "Enter quarter number to fill: \n";
                cin>> quarter;
                fillingWithLines(hdc, xc1, yc1, r1, color(colorChoice),quarter);
                quarter=0;
                break;
            case FillCircle_wcircle:
                cout<< "Enter quarter number to fill: \n";
                cin>> quarter;
                fillingWithCircles(hdc, xc1, yc1, r1, color(colorChoice),quarter);
                quarter=0;
                break;
            // recursive flood fill
            case FillRecursive:
                choice = FillRecursive;
                break ;
            // non_recursive flood fill
            case FillnonRecursive:
                choice = FillnonRecursive;
                break ;
            //square filling
             case FillSquare:
                choice = square ;
                break;
            //rectangle filling
            case FillRectangle:
                choice = rectangle;
                break;
            //convex filling
            case FillConvex:
                choice = FillConvex;
                copyPoints();
                ConvexFill(hdc, pointArr , counter , color(colorChoice));
                break ;
            //non convex filling
            case FillnonConvex :
                choice = FillnonConvex;
                copyPoints();
                GeneralPolygonFill(hdc,pointArr,counter,color(colorChoice));
                break ;
            //recursive filling

            //non recursive filling

            }
        break ;
        /*clicks*/
        case WM_LBUTTONDOWN:
            if (choice >= 9&&choice<=11){
                    if (t == 1){
                    x1 = LOWORD(lParam);
                    y1f = HIWORD(lParam);
                    dots[counter].x = x1 ;
                    dots[counter].y = y1f;
                    counter ++ ;
                }else if(t==2){
                    x2 = LOWORD(lParam);
                    y2 = HIWORD(lParam);
                    dots[counter].x = x2 ;
                    dots[counter].y = y2 ;
                    counter ++ ;
                    t = 0;
                    if(choice==9){
                        DrawLineDDA(hdc, x1, y1f, x2, y2, color(colorChoice));
                    }else if(choice==10){
                        DrawLineMidpoint(hdc, x1, y1f, x2, y2, color(colorChoice));
                    }else if(choice==11){
                        DrawLineParametric(hdc, x1, y1f, x2, y2, color(colorChoice));
                    }
                }
            }else if(choice>=12&&choice<=16){
                if (t == 1){
                    xc1 = LOWORD(lParam);
                    yc1 = HIWORD(lParam);
                }else if(t==2){
                    xc2 = LOWORD(lParam);
                    yc2 = HIWORD(lParam);
                    t = 0;
                    r1 = sqrt(pow(xc1-xc2,2.0)+pow(yc1-yc2,2.0));
                    if(choice==12){
                        DrawCircleDirect(hdc, xc1, yc1, r1, color(colorChoice),quarter);
                        quarter=0;
                    }else if(choice==13){
                        DrawCirclePolar(hdc, xc1, yc1, r1, color(colorChoice),quarter);
                        quarter=0;
                    }else if(choice==14){
                        DrawCircleIterativePolar(hdc, xc1, yc1, r1, color(colorChoice),quarter);
                        quarter=0;
                    }else if(choice==15){
                        DrawCircleMidPoint(hdc, xc1, yc1, r1, color(colorChoice),quarter);
                        quarter=0;
                    }else if(choice==16){
                        DrawCircleMidPoint_modified(hdc, xc1, yc1, r1, color(colorChoice),quarter);
                        quarter=0;
                    }
                }
            }else if (choice >= 17&&choice<=19){
                if (t == 1){
                    x1 = LOWORD(lParam);
                    y1f = HIWORD(lParam);
                }else if (t==2){
                    x2 = LOWORD(lParam);
                    y2 = HIWORD(lParam);
                    A = calcR(x1,y1f,x2,y2);
                }else if(t==3){
                    x3 = LOWORD(lParam);
                    y3 = HIWORD(lParam);
                    B = calcR(x1,y1f,x3,y3);
                    t = 0;
                    if(choice==17){
                        DrawEllipseDirect(hdc, x1, y1f, A, B, color(colorChoice));
                    }else if(choice==18){
                        DrawEllipsePolar(hdc, x1, y1f, A, B, color(colorChoice));
                    }else if(choice==19){
                        midellipse(hdc, x1, y1f, A, B, color(colorChoice));
                    }
                }
            }else if (choice >= 20&&choice<=26){
                if(choice >= 20&&choice<=24){
                    if (t == 1){
                        X_left = LOWORD(lParam);
                        Y_top = HIWORD(lParam);
                    }else if (t == 2){
                        X_right = LOWORD(lParam);
                        Y_bottom = HIWORD(lParam);
                        Rectangle(hdc, X_left, Y_top, X_right, Y_bottom);
                    }else if(t == 3){
                        x1 = LOWORD(lParam);
                        yy1 = HIWORD(lParam);
                        if(choice==20||choice==23){
                            PointClipping(hdc, x1, yy1, X_left, Y_top, X_right, Y_bottom);
                            t--;
                        }
                    }else if (t == 4) {
                        x2 = LOWORD(lParam);
                        y2 = HIWORD(lParam);
                        if(choice==21||choice==24){
                            CohenSuth(hdc, x1, yy1, x2, y2, X_left, Y_top, X_right, Y_bottom);
                            t -=2;
                        }
                    }else if(choice==22){
                        int n = 5;
                        int x = LOWORD(lParam);
                        int y = HIWORD(lParam);
                        p1[t - 2] = Vertex((double)x, (double)y);

                        if (t == 7) {
                            PolygonClip(hdc, p1, n, X_left, Y_top, X_right, Y_bottom);
                            t = 2;
                        }

                    }

                }else if(choice == 25||choice==26){
                    if (t == 1){
                        xc1 = LOWORD(lParam);
                        yc1 = HIWORD(lParam);
                    }else if(t==2){
                        xc2 = LOWORD(lParam);
                        yc2 = HIWORD(lParam);
                        r1 = sqrt(pow(xc1-xc2,2.0)+pow(yc1-yc2,2.0));
                        data.push_back({22,xc1, yc1, r1, colorChoice,quarter});
                        DrawCircleDirect(hdc, xc1, yc1, r1, color(colorChoice),quarter);
                        cout<<t;
                    }else if(t == 3) {
                        x1 = LOWORD(lParam);
                        yy1 = HIWORD(lParam);
                        if(choice==25){
                            data.push_back({25, x1, yy1, xc1, yc1, r1});
                            PointClippingInCircle(hdc, x1, yy1, xc1, yc1, r1);
                            cout<<t;
                            t--;
                        }
                    }else if (t == 4) {
                        x2 = LOWORD(lParam);
                        y2 = HIWORD(lParam);
                        if(choice==26){
                            data.push_back({26,x1, x2, yy1, y2, xc1, yc1, r1 });
                            LineClippingInCircle(hdc, x1, x2, yy1, y2, xc1, yc1, r1);
                            cout<<t;
                            t-=2;
                        }
                    }
                }

            }else if(choice == 33||choice == 34){
                x1 = LOWORD(lParam);
                y1f = HIWORD(lParam);
                if(choice == 33){
                    FloodFill(hdc , x1,y1f,color(colorChoice),color(colorChoice));
                }else if(choice == 34){
                    NRFloodFill(hdc , x1,y1f,color(colorChoice),color(colorChoice));
                }
                t--;
            }else if (choice == FillConvex){
                t--;
            }
            else if (choice == FillnonConvex){
                t-- ;
            }
            else if (choice == square){
                if (t == 1){
                    x1 = LOWORD(lParam);
                    y1f = HIWORD(lParam);
                }else if (t == 2){
                    x2 = LOWORD(lParam);
                    y2 = HIWORD(lParam);
                    Rectangle(hdc, x1, y1f, x2, y2);
                }else if ( t == 3 ){
                    t = 0 ;
                    for (int i = x1; i <= x2; i++)
                    {
                    p[0] = Vector(i, y1f);
                    p[1] = Vector(i, y1f - 1);
                    p[2] = Vector(i, y1f - 2);
                    p[3] = Vector(i, y2);
                    Vector T1(3 * (p[1][0] - p[0][0]), 3 * (p[1][1] - p[0][1]));
                    Vector T2(3 * (p[3][0] - p[2][0]), 3 * (p[3][1] - p[2][1]));
                    DrawHermiteCurve(hdc, p[0], T1, p[3], T2,color(colorChoice));
                    }
                }

            }else if (choice == rectangle){
                if (t == 1){
                    x1 = LOWORD(lParam);
                    y1f = HIWORD(lParam);
                }else if (t == 2){
                    x2 = LOWORD(lParam);
                    y2 = HIWORD(lParam);
                    Rectangle(hdc, x1, y1f, x2, y2);
                }else if (t == 3){
                    t = 0 ;
                    for (int i = y1f; i <= y2; i++){
                    p[0] = Vector(x1, i);
                    p[1] = Vector(x1 - 1,i );
                    p[2] = Vector( x1 - 2,i);
                    p[3] = Vector(x2, i);
                    DrawBezierCurve(hdc, p[0], p[1], p[2], p[3],color(colorChoice));
                    }
                }
            }else if (choice == SplineCurve){
                if(t == 6){
                    t=0;
                    DrawCardinalSpline(hdc,splinepoint,6,0.2,color(colorChoice));
                }else{
                    splinepoint[t].v[0] = LOWORD(lParam);
                    splinepoint[t].v[1] = HIWORD(lParam);
                }
            }
            t++;
            break;
        case WM_RBUTTONDOWN :
            counter = 0 ;
            break ;
        case WM_CREATE:
            hdc = GetDC(hwnd);
            AddMenus(hwnd);
            break;
        case WM_DESTROY:
            PostQuitMessage (0);       /* send a WM_QUIT to the message queue */
            break;
        default:                      /* for messages that we don't deal with */
        return DefWindowProc (hwnd, message, wParam, lParam);
    }

    return 0;
}
void BKcolor(HWND hwnd,COLORREF BKc){
    HBRUSH hBrush = CreateSolidBrush(BKc);
    SetClassLongPtr(hwnd, GCLP_HBRBACKGROUND, HandleToLong(hBrush));
    InvalidateRect(hwnd, NULL, TRUE);
}

void AddMenus(HWND hwnd){
    hmenu = CreateMenu();

    HMENU hFileMenu = CreateMenu();
    HMENU hBKMenu = CreateMenu();

    HMENU hDrawMenu = CreateMenu();
    HMENU hDrawLineMenu = CreateMenu();
    HMENU hDrawCircleMenu = CreateMenu();
    HMENU hDrawEllipseMenu = CreateMenu();

    HMENU hClippingMenu = CreateMenu();
    HMENU hClipp_Rectangle = CreateMenu();
    HMENU hClipp_Square = CreateMenu();
    HMENU hClipp_Circle = CreateMenu();


    HMENU hFillMenu = CreateMenu();
    HMENU hFillCircleMenu = CreateMenu();

    AppendMenu(hmenu,MF_POPUP,(UINT_PTR)hFileMenu,"File");
    AppendMenu(hFileMenu,MF_STRING,File_save,"Save");
    AppendMenu(hFileMenu,MF_STRING,File_load,"Load");
    AppendMenu(hFileMenu,MF_STRING,File_choseColor,"Chose Color");
    AppendMenu(hFileMenu,MF_POPUP,(UINT_PTR)hBKMenu,"Background Color");
        AppendMenu(hBKMenu,MF_STRING,BK_color_white,"Background Color White");
        AppendMenu(hBKMenu,MF_STRING,BK_color_gray,"Background Color Gray");
        AppendMenu(hBKMenu,MF_STRING,BK_color_black,"Background Color Black");
    AppendMenu(hFileMenu,MF_SEPARATOR,NULL,NULL);
    AppendMenu(hFileMenu,MF_STRING,File_clear,"Clear");
    AppendMenu(hFileMenu,MF_STRING,File_exit,"Exit");

    AppendMenu(hmenu,MF_POPUP,(UINT_PTR)hDrawMenu,"Draw");

    AppendMenu(hDrawMenu,MF_POPUP,(UINT_PTR)hDrawLineMenu,"Line");
    AppendMenu(hDrawLineMenu,MF_STRING,DrawLine_dda,"DDA");
    AppendMenu(hDrawLineMenu,MF_STRING,DrawLine_midline,"Midpoint");
    AppendMenu(hDrawLineMenu,MF_STRING,DrawLine_paramline,"Parametric");

    AppendMenu(hDrawMenu,MF_POPUP,(UINT_PTR)hDrawCircleMenu,"Circle");
    AppendMenu(hDrawCircleMenu,MF_STRING,DrawCircle_Direct,"Direct");
    AppendMenu(hDrawCircleMenu,MF_STRING,DrawCircle_Polar,"Polar");
    AppendMenu(hDrawCircleMenu,MF_STRING,DrawCircle_ItPolar,"Iterative Polar");
    AppendMenu(hDrawCircleMenu,MF_STRING,DrawCircle_midcircle,"Midpoint");
    AppendMenu(hDrawCircleMenu,MF_STRING,DrawCircle_midcircleModified,"Modified Midpoint");

    AppendMenu(hDrawMenu,MF_POPUP,(UINT_PTR)hDrawEllipseMenu,"Ellipse");
    AppendMenu(hDrawEllipseMenu,MF_STRING,DrawEllipse_Direct,"Direct");
    AppendMenu(hDrawEllipseMenu,MF_STRING,DrawEllipse_Polar,"Polar");
    AppendMenu(hDrawEllipseMenu,MF_STRING,DrawEllipse_midellipse,"Midpoint");

    AppendMenu(hmenu,MF_POPUP,(UINT_PTR)hClippingMenu,"Clipping");
    AppendMenu(hClippingMenu,MF_POPUP,(UINT_PTR)hClipp_Rectangle,"Rectangle");
        AppendMenu(hClipp_Rectangle,MF_STRING,Clipp_RecPoint,"Point");
        AppendMenu(hClipp_Rectangle,MF_STRING,Clipp_RecLine,"Line");
        AppendMenu(hClipp_Rectangle,MF_STRING,Clipp_RecPolygon,"Polygon");
    AppendMenu(hClippingMenu,MF_POPUP,(UINT_PTR)hClipp_Square,"Square");
        AppendMenu(hClipp_Square,MF_STRING,Clipp_SquPoint,"Point");
        AppendMenu(hClipp_Square,MF_STRING,Clipp_SquLine,"Line");
    AppendMenu(hClippingMenu,MF_POPUP,(UINT_PTR)hClipp_Circle,"Circle");
        AppendMenu(hClipp_Circle,MF_STRING,Clipp_CirclePoint,"Point");
        AppendMenu(hClipp_Circle,MF_STRING,Clipp_CircleLine,"Line");

    AppendMenu(hmenu,MF_POPUP,(UINT_PTR)hFillMenu,"Fill");
    AppendMenu(hFillMenu,MF_POPUP,(UINT_PTR)hFillCircleMenu,"Circle");
        AppendMenu(hFillCircleMenu,MF_STRING,FillCircle_wline,"With Lines");
        AppendMenu(hFillCircleMenu,MF_STRING,FillCircle_wcircle,"With Circles");
    AppendMenu(hFillMenu,MF_STRING,FillSquare,"Square");
    AppendMenu(hFillMenu,MF_STRING,FillRectangle,"Rectangle");
    AppendMenu(hFillMenu,MF_SEPARATOR,NULL,NULL);
    AppendMenu(hFillMenu,MF_STRING,FillConvex,"Convex");
    AppendMenu(hFillMenu,MF_STRING,FillnonConvex,"Non Convex");
    AppendMenu(hFillMenu,MF_SEPARATOR,NULL,NULL);
    AppendMenu(hFillMenu,MF_STRING,FillRecursive,"Recursive");
    AppendMenu(hFillMenu,MF_STRING,FillnonRecursive,"Non Recursive");

    AppendMenu(hmenu,MF_STRING,SplineCurve,"Spline Curve");
    SetMenu(hwnd,hmenu);
}

