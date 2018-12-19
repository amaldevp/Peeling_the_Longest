#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include<string.h>
#include<fstream>
#include<string.h>
#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#define NENDS 2
float new_pts[100000][3];
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3   Point;
GLdouble width, height;
int wd,loop,ends1[NENDS][2],di=0;
Delaunay dt;
float len;
Point deleted_pts[100000],saved_pts[10000],input_points[10000];
int dpi=0,spi=0;
int vertex_count=0,minx=9999,miny=9999,maxx=0,maxy=0,vertex_size;
void init(void)
{
    width  = 1280.0;
    height = 800.0;
    ends1[0][0] = (int)(0.25*width);
    ends1[0][1] = (int)(0.75*height);
    ends1[1][0] = (int)(0.75*width);
    ends1[1][1] = (int)(0.25*height);
    return;
}
void reshape(int w, int h)
{
    width = (GLdouble) w;
    height = (GLdouble) h;
    glViewport(0, 0, (GLsizei) width, (GLsizei) height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glRotatef(-180,180,0,1);
    glOrtho(minx-20.0,maxx+20.0,miny-20.0,maxy+20.0, -1.f, 1.f);
    return;
}
void kbd(unsigned char key, int x, int y)
{
    switch((char)key) {
        case 'q':
        case 27:
            glutDestroyWindow(wd);
            exit(0);
        default:
            break;
    }
    return;
}
void drawFilledCircle(GLfloat x, GLfloat y, GLfloat radius)
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    int i,triangleAmount = 20;
    GLfloat twicePi = 2.0f * 3.14159265;
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y);
    for(i = 0; i <= triangleAmount;i++)
        glVertex2f(x + (radius * cos(i *  twicePi / triangleAmount)),y + (radius * sin(i * twicePi / triangleAmount)));
    glEnd();
}
char *fname;
void pointset(void)
{
    glLineWidth(6.0);
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_POLYGON_SMOOTH );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
    glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );
    glColor3f(0,0,0);
    for(int i=0;i<spi;i++)
        drawFilledCircle(saved_pts[i].x(),saved_pts[i].y(),vertex_size);
}
void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);
    pointset();
    glFlush();
    return;
}
float distance(Point a, Point b)
{
    return (float)(sqrt(abs(((a.x()-b.x())*(a.x()-b.x())))+abs(((a.y()-b.y())*(a.y()-b.y())))));
}
void delete_pts(Point a)
{
    saved_pts[spi]=a;
    spi++;
    for (Delaunay::Finite_vertices_iterator it = dt.finite_vertices_begin(); it != dt.finite_vertices_end(); it++)
        if(distance(it->point(),a)<len)
        {
            deleted_pts[dpi]=it->point();
            dpi++;
        }
}
int check_in(Point a)
{
    for(int i=0;i<dpi;i++)
        if(a==deleted_pts[i])
            return 1;
    return 0;
}
int main(int argc, char **argv)
{
    init();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowSize(250, 250);
    wd = glutCreateWindow("Peel the longest: Noise removal");
    repeat: std::ifstream ifs(argv[1]);
    if(ifs.fail())
    {
        printf("File does not exists\n");
        exit(1);
    }
    std::istream_iterator<Point> begin(ifs);
    std::istream_iterator<Point> end;
    dt.insert(begin, end);
    Delaunay::Vertex_iterator vi=dt.vertices_begin();
    ifs.close();
    do{
        if(vi->point().x()>maxx)
            maxx=vi->point().x();
        if(vi->point().y()>maxy)
            maxy=vi->point().y();
        if(vi->point().x()<minx)
            minx=vi->point().x();
        if(vi->point().y()<miny)
            miny=vi->point().y();
        input_points[vertex_count]=vi->point();
        vertex_count++;
        if(vertex_count>10000)
        {
            printf("Input size is too high, please increase the size of input_points and shape array\n");
            exit(1);
        }
        vi++;
    }while(vi!=dt.vertices_end());
    printf("Enter value for noise tuning parameter: ");
    scanf("%f",&len);
    for (Delaunay::Finite_vertices_iterator it = dt.finite_vertices_begin(); it != dt.finite_vertices_end(); it++)
    {
        if(check_in(it->point()))
            continue;
        float dist[10000];
        int di=0;
        Delaunay::Vertex_circulator circulator = dt.incident_vertices(it), done(circulator);
        do
        {
            dist[di]=distance(circulator->point(),it->point());
            di++;
        } while(++circulator != done);
        float big=0,small=9999;
        for(int i=0;i<di;i++)
        {
            if(dist[i]>big)
            {
                big=dist[i];
            }
            if(dist[i]<small)
                small=dist[i];
        }
        if(big<2*small)
            delete_pts(it->point());
    }
    printf("Enter vertex size for displaying: ");
    scanf("%d",&vertex_size);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(kbd);
    glutDisplayFunc(display);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(3.0);
    glutMainLoop();
    exit(0);
    return 0;
}
