#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include<string.h>
#include<fstream>
#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#define NENDS 2
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3   Point;
GLdouble width, height;
int wd,loop,ends1[NENDS][2],di=0;
Point shape[500000][2],input_points[500000],deleted[500000][2],shape1[500000][2],shape2[500000][2],tempv[1000][2];
int vertex_count=0,shape_index=0,minx=9999,miny=9999,maxx=0,maxy=0,vertex_size,shape_fl,shape1_index=0,shape2_index=0;
Delaunay dt;
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
void pointset(void)
{
    glLineWidth(6.0);
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_POLYGON_SMOOTH );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
    glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );
    for(int i=0;i<shape1_index;i++)
    {
        glColor3f(0.5,0.5,0.5);
        glBegin(GL_LINES);
        glVertex2f(shape1[i][0].x(),shape1[i][0].y());
        glVertex2f(shape1[i][1].x(),shape1[i][1].y());
        glEnd();
    }
    glColor3f(0,0,0);
    for(int i=0;i<vertex_count;i++)
        drawFilledCircle(input_points[i].x(),input_points[i].y(),vertex_size);
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
void insert_to_shape(Point a,Point b)
{
    for(int i=0;i<shape_index;i++)
        if((shape[i][0]==a&&shape[i][1]==b)||(shape[i][0]==b&&shape[i][1]==a))
            return;
    shape[shape_index][0]=a;
    shape[shape_index][1]=b;
    shape_index++;
}
bool notsmallest(Delaunay::Vertex_handle a,Delaunay::Vertex_handle b)
{
    Delaunay::Vertex_handle vh1;
    Delaunay::Vertex_circulator vc=dt.incident_vertices(a),done(vc);
    double min1=99999;
    if (vc != 0) {
        do {
            if(distance(vc->point(),a->point())<min1)
            {
                min1=distance(vc->point(),a->point());
                vh1=vc;
            }
        }while(++vc != done);
        if(vh1->point()==b->point())
            return 0;
    }
    Delaunay::Vertex_circulator vc1=dt.incident_vertices(b),done1(vc1);
    min1=99999;
    if (vc1 != 0) {
        do {
            if(distance(vc1->point(),b->point())<min1)
            {
                min1=distance(vc1->point(),b->point());
                vh1=vc1;
            }
        }while(++vc1 != done1);
        if(vh1->point()==a->point())
            return 0;
    }
    return 1;
}
int search_in_shape(Point a,Point b)
{
    for(int i=0;i<shape_index;i++)
        if((shape[i][0]==a&&shape[i][1]==b)||(shape[i][0]==b&&shape[i][1]==a))
            return i;
    return -1;
}
int deleted_edge(Point a,Point b)
{
    for(int i=0;i<di;i++)
        if((deleted[i][0]==a&&deleted[i][1]==b)||(deleted[i][0]==b&&deleted[i][1]==a))
            return 1;
    return 0;
}
int notexist(Point a,Point b)
{
    for(int i=0;i<shape1_index;i++)
    {
        if(shape1[i][0]==a&&shape1[i][1]==b)
            return 0;
        if(shape1[i][0]==b&&shape1[i][1]==a)
            return 0;
    }
    return 1;
}
int main(int argc, char **argv)
{
    init();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowSize(250, 250);
    wd = glutCreateWindow("Crawl through Neighbors");
    std::ifstream ifs(argv[1]);
    if(ifs.fail())
    {
        printf("File does not exists\n");
        exit(1);
    }
    FILE *nfp1;
    nfp1=fopen("new.txt","w");
    std::istream_iterator<Point> begin(ifs);
    std::istream_iterator<Point> end;
    dt.insert(begin, end);
    Delaunay::Vertex_iterator vi=dt.vertices_begin();
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
        fprintf(nfp1,"%d,%d\n",(int)vi->point().x(),(int)vi->point().y());
        if(vertex_count>500000)
        {
            printf("Input size is too high, please increase the size of input_points and shape array\n");
            exit(1);
        }
        vi++;
    }while(vi!=dt.vertices_end());
    fclose(nfp1);
    printf("Enter the parameter for restoring self-intersections (0 if no self-intersections): ");
    float len;
    scanf("%f",&len);
    for (Delaunay::Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++)
    {
        if(distance(dt.triangle(it)[0],dt.triangle(it)[1])>=distance(dt.triangle(it)[0],dt.triangle(it)[2])
           &&distance(dt.triangle(it)[0],dt.triangle(it)[1])>=distance(dt.triangle(it)[1],dt.triangle(it)[2]))
        {
            deleted[di][0]=dt.triangle(it)[0];
            deleted[di][1]=dt.triangle(it)[1];
            di++;
            int index=search_in_shape(dt.triangle(it)[0],dt.triangle(it)[1]);
            if(index>-1)
            {
                for(int i=index;i<shape_index-1;i++)
                {
                    shape[i][0]=shape[i+1][0];
                    shape[i][1]=shape[i+1][1];
                }
                shape_index--;
            }
            if(!deleted_edge(dt.triangle(it)[0],dt.triangle(it)[2]))
                insert_to_shape(dt.triangle(it)[0],dt.triangle(it)[2]);
            if(!deleted_edge(dt.triangle(it)[1],dt.triangle(it)[2]))
                insert_to_shape(dt.triangle(it)[2],dt.triangle(it)[1]);
        }
        if(distance(dt.triangle(it)[0],dt.triangle(it)[2])>=distance(dt.triangle(it)[0],dt.triangle(it)[1])
           &&distance(dt.triangle(it)[0],dt.triangle(it)[2])>=distance(dt.triangle(it)[1],dt.triangle(it)[2]))
        {
            deleted[di][0]=dt.triangle(it)[0];
            deleted[di][1]=dt.triangle(it)[2];
            di++;
            int index=search_in_shape(dt.triangle(it)[0],dt.triangle(it)[2]);
            if(index>-1)
            {
                for(int i=index;i<shape_index-1;i++)
                {
                    shape[i][0]=shape[i+1][0];
                    shape[i][1]=shape[i+1][1];
                }
                shape_index--;
            }
            if(!deleted_edge(dt.triangle(it)[0],dt.triangle(it)[1]))
                insert_to_shape(dt.triangle(it)[0],dt.triangle(it)[1]);
            if(!deleted_edge(dt.triangle(it)[1],dt.triangle(it)[2]))
                insert_to_shape(dt.triangle(it)[1],dt.triangle(it)[2]);
        }
        if(distance(dt.triangle(it)[2],dt.triangle(it)[1])>=distance(dt.triangle(it)[0],dt.triangle(it)[1])
           &&distance(dt.triangle(it)[2],dt.triangle(it)[1])>=distance(dt.triangle(it)[0],dt.triangle(it)[2]))
        {
            deleted[di][0]=dt.triangle(it)[2];
            deleted[di][1]=dt.triangle(it)[1];
            di++;
            int index=search_in_shape(dt.triangle(it)[2],dt.triangle(it)[1]);
            if(index>-1)
            {
                for(int i=index;i<shape_index-1;i++)
                {
                    shape[i][0]=shape[i+1][0];
                    shape[i][1]=shape[i+1][1];
                }
                shape_index--;
            }
            if(!deleted_edge(dt.triangle(it)[0],dt.triangle(it)[1]))
                insert_to_shape(dt.triangle(it)[0],dt.triangle(it)[1]);
            if(!deleted_edge(dt.triangle(it)[0],dt.triangle(it)[2]))
                insert_to_shape(dt.triangle(it)[0],dt.triangle(it)[2]);
        }
    }
    for(int i=0;i<shape_index;i++)
    {
        Point first[1000];
        int fi=0;
        float first_dist[1000];
        for(int j=0;j<shape_index;j++)
        {
            if(j!=i)
                if(shape[i][0]==shape[j][0]||shape[i][0]==shape[j][1])
                {
                    if(shape[i][0]==shape[j][0])
                        first[fi]=shape[j][1];
                    if(shape[i][0]==shape[j][1])
                        first[fi]=shape[j][0];
                    fi++;
                }
        }
        int degree_count=0;
        float smallest=9999,smallest1=9999;
        if(fi>1)
        {
            for(int j=0;j<fi;j++)
                first_dist[j]=distance(shape[i][0],first[j]);
            float org_dist=distance(shape[i][0],shape[i][1]);
            smallest=first_dist[0];
            for(int j=0;j<fi;j++)
                if(first_dist[j]<org_dist)
                {
                    if(first_dist[j]<smallest)
                        smallest=first_dist[j];
                    degree_count++;
                }
        }
        for(int j=0;j<shape_index;j++)
        {
            if(j!=i)
                if(shape[i][1]==shape[j][0]||shape[i][1]==shape[j][1])
                {
                    if(shape[i][1]==shape[j][0])
                        first[fi]=shape[j][1];
                    if(shape[i][1]==shape[j][1])
                        first[fi]=shape[j][0];
                    fi++;
                }
        }
        int degree_count1=0;
        if(fi>1)
        {
            for(int j=0;j<fi;j++)
                first_dist[j]=distance(shape[i][1],first[j]);
            float org_dist=distance(shape[i][0],shape[i][1]);
            smallest1=first_dist[0];
            for(int j=0;j<fi;j++)
                if(first_dist[j]<org_dist)
                {
                    if(first_dist[j]<smallest1)
                        smallest1=first_dist[j];
                    degree_count1++;
                }
        }
        if(smallest1>smallest)
            smallest=smallest1;
        if((degree_count<2&&degree_count1<2))
        {
            shape1[shape1_index][0]=shape[i][0];
            shape1[shape1_index][1]=shape[i][1];
            shape1_index++;
        }
    }
    for(int i=0;i<shape1_index;i++)
    {
        int deg1=0,deg2=0;
        for(int j=0;j<shape1_index;j++)
        {
            if(j!=i)
            {
                if(shape1[i][1]==shape1[j][0]||shape1[i][1]==shape1[j][1])
                {
                    deg1++;
                }
                if(shape1[i][0]==shape1[j][0]||shape1[i][0]==shape1[j][1])
                {
                    deg2++;
                }
            }
        }
        Point tv1,tv2;
        double min=9999;
        int fl=0;
        if(((deg2==0)))
            for (Delaunay::Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++)
            {
                if(dt.triangle(it)[0]==shape1[i][0])
                {
                    if(distance(dt.triangle(it)[0],dt.triangle(it)[1])<len*distance(shape1[i][0],shape1[i][1]))
                    {
                        if(min>distance(dt.triangle(it)[0],dt.triangle(it)[1])&&notexist(dt.triangle(it)[0],dt.triangle(it)[1]))
                        {
                            min=distance(dt.triangle(it)[0],dt.triangle(it)[1]);
                            tv1=dt.triangle(it)[1];
                            tv2=dt.triangle(it)[0];
                            fl=1;
                        }
                        
                    }
                    if(distance(dt.triangle(it)[0],dt.triangle(it)[2])<len*distance(shape1[i][0],shape1[i][1]))
                    {
                        if(min>distance(dt.triangle(it)[0],dt.triangle(it)[2])&&notexist(dt.triangle(it)[0],dt.triangle(it)[2]))
                        {
                            min=distance(dt.triangle(it)[0],dt.triangle(it)[2]);
                            tv1=dt.triangle(it)[2];
                            tv2=dt.triangle(it)[0];
                            fl=1;
                        }
                    }
                }
                else
                {
                    if(dt.triangle(it)[1]==shape1[i][0])
                    {
                        if(distance(dt.triangle(it)[0],dt.triangle(it)[1])<len*distance(shape1[i][0],shape1[i][1]))
                        {
                            if(min>distance(dt.triangle(it)[0],dt.triangle(it)[1])&&notexist(dt.triangle(it)[0],dt.triangle(it)[1]))
                            {
                                min=distance(dt.triangle(it)[0],dt.triangle(it)[1]);
                                tv1=dt.triangle(it)[1];
                                tv2=dt.triangle(it)[0];
                                fl=1;
                            }
                        }
                        if(distance(dt.triangle(it)[2],dt.triangle(it)[1])<len*distance(shape1[i][0],shape1[i][1]))
                        {
                            if(min>distance(dt.triangle(it)[2],dt.triangle(it)[1])&&notexist(dt.triangle(it)[2],dt.triangle(it)[1]))
                            {
                                min=distance(dt.triangle(it)[2],dt.triangle(it)[1]);
                                tv1=dt.triangle(it)[1];
                                tv2=dt.triangle(it)[2];
                                fl=1;
                            }
                        }
                    }
                    else
                    {
                        if(dt.triangle(it)[2]==shape1[i][0])
                        {
                            if(distance(dt.triangle(it)[0],dt.triangle(it)[2])<len*distance(shape1[i][0],shape1[i][1]))
                            {
                                if(min>distance(dt.triangle(it)[0],dt.triangle(it)[2])&&notexist(dt.triangle(it)[0],dt.triangle(it)[2]))
                                {
                                    min=distance(dt.triangle(it)[0],dt.triangle(it)[2]);
                                    tv1=dt.triangle(it)[2];
                                    tv2=dt.triangle(it)[0];
                                    fl=1;
                                }
                            }
                            if(distance(dt.triangle(it)[2],dt.triangle(it)[1])<len*distance(shape1[i][0],shape1[i][1]))
                            {
                                if(min>distance(dt.triangle(it)[1],dt.triangle(it)[2])&&notexist(dt.triangle(it)[2],dt.triangle(it)[1]))
                                {
                                    min=distance(dt.triangle(it)[1],dt.triangle(it)[2]);
                                    tv1=dt.triangle(it)[1];
                                    tv2=dt.triangle(it)[2];
                                    fl=1;
                                }
                                
                            }
                        }
                        
                    }
                }
                
            }
        if(fl==1)
        {
            shape1[shape1_index][0]=tv1;
            shape1[shape1_index][1]=tv2;
            shape1_index++;
        }
        min=9999;
        fl=0;
        if(((deg1==0)))
            for (Delaunay::Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++)
            {
                if(dt.triangle(it)[0]==shape1[i][1])
                {
                    if(distance(dt.triangle(it)[0],dt.triangle(it)[1])<len*distance(shape1[i][0],shape1[i][1]))
                    {
                        if(min>distance(dt.triangle(it)[0],dt.triangle(it)[1])&&notexist(dt.triangle(it)[0],dt.triangle(it)[1]))
                        {
                            min=distance(dt.triangle(it)[0],dt.triangle(it)[1]);
                            tv1=dt.triangle(it)[1];
                            tv2=dt.triangle(it)[0];
                            fl=1;
                        }
                        
                    }
                    if(distance(dt.triangle(it)[0],dt.triangle(it)[2])<len*distance(shape1[i][0],shape1[i][1]))
                    {
                        if(min>distance(dt.triangle(it)[0],dt.triangle(it)[2])&&notexist(dt.triangle(it)[0],dt.triangle(it)[2]))
                        {
                            min=distance(dt.triangle(it)[0],dt.triangle(it)[2]);
                            tv1=dt.triangle(it)[2];
                            tv2=dt.triangle(it)[0];
                            fl=1;
                        }
                    }
                }
                else
                {
                    if(dt.triangle(it)[1]==shape1[i][1])
                    {
                        if(distance(dt.triangle(it)[0],dt.triangle(it)[1])<len*distance(shape1[i][0],shape1[i][1]))
                        {
                            if(min>distance(dt.triangle(it)[0],dt.triangle(it)[1])&&notexist(dt.triangle(it)[0],dt.triangle(it)[1]))
                            {
                                min=distance(dt.triangle(it)[0],dt.triangle(it)[1]);
                                tv1=dt.triangle(it)[1];
                                tv2=dt.triangle(it)[0];
                                fl=1;
                            }
                        }
                        if(distance(dt.triangle(it)[2],dt.triangle(it)[1])<len*distance(shape1[i][0],shape1[i][1]))
                        {
                            if(min>distance(dt.triangle(it)[2],dt.triangle(it)[1])&&notexist(dt.triangle(it)[2],dt.triangle(it)[1]))
                            {
                                min=distance(dt.triangle(it)[2],dt.triangle(it)[1]);
                                tv1=dt.triangle(it)[1];
                                tv2=dt.triangle(it)[2];
                                fl=1;
                            }
                        }
                    }
                    else
                    {
                        if(dt.triangle(it)[2]==shape1[i][1])
                        {
                            if(distance(dt.triangle(it)[0],dt.triangle(it)[2])<len*distance(shape1[i][0],shape1[i][1]))
                            {
                                if(min>distance(dt.triangle(it)[0],dt.triangle(it)[2])&&notexist(dt.triangle(it)[0],dt.triangle(it)[2]))
                                {
                                    min=distance(dt.triangle(it)[0],dt.triangle(it)[2]);
                                    tv1=dt.triangle(it)[2];
                                    tv2=dt.triangle(it)[0];
                                    fl=1;
                                }
                            }
                            if(distance(dt.triangle(it)[2],dt.triangle(it)[1])<len*distance(shape1[i][0],shape1[i][1]))
                            {
                                if(min>distance(dt.triangle(it)[1],dt.triangle(it)[2])&&notexist(dt.triangle(it)[2],dt.triangle(it)[1]))
                                {
                                    min=distance(dt.triangle(it)[1],dt.triangle(it)[2]);
                                    tv1=dt.triangle(it)[1];
                                    tv2=dt.triangle(it)[2];
                                    fl=1;
                                }
                                
                            }
                        }
                        
                    }
                }
                
            }
        if(fl==1)
        {
            shape1[shape1_index][0]=tv1;
            shape1[shape1_index][1]=tv2;
            shape1_index++;
        }
    }
    printf("Enter vertex size\n");
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
