#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <conio.h>

using namespace std;

struct Joint
{
    double x;
    double y;
    int ID;
};

int Original_Truss_Color;
string Analysed_Truss_Show_Option;
double Members_fy;

struct Member
{
    int i;
    int j;
    int ID;
    double l;
    double l2;
    double Sin;
    double Cos;
    double A;
    double E;
    double AE;
    double Stress;
    double Delta;
    int Color;

};

struct Support
{
    int ID;
    int j;
    int xR;
    int yR;
};

struct Load
{
    int ID;
    int j;
    double xC;
    double yC;
    double Mag;
};

void MakeDXF( Joint *Joints, Joint *FinalJoints, Member *Members, Support *Supports,Load *Loads , int MemberColor, int SupportColor, int FinalSupportColor, int LoadColor, int JointCount, int MemberCount, int SupportCount, int LoadCount );

int main()
{

    // Declare File Object
    fstream InputFile;

    // Open Input File
    InputFile.open("TA.input", ios::in);

    //Input Original Truss Color
    InputFile >> Original_Truss_Color;

    //Input Analysed Truss Show Option
    InputFile >> Analysed_Truss_Show_Option;

    //Input Joints
    int JointNumbers;
    InputFile >> JointNumbers;
    Joint Joints[JointNumbers];
    for (int i = 0; i < JointNumbers; ++i)
    {
        Joints[i].ID = i + 1;
        InputFile >> Joints[i].x;
        InputFile >> Joints[i].y;
    }

    //Input Members
    int MemberNumbers;
    InputFile >> MemberNumbers;
    InputFile >> Members_fy;
    Member Members[MemberNumbers];
    for (int i = 0; i < MemberNumbers ; ++i)
    {
        Members[i].ID = i + 1;
        InputFile >> Members[i].i;
        InputFile >> Members[i].j;
        InputFile >> Members[i].A;
        InputFile >> Members[i].E;
    }

    // Input Supports
    int SupportNumbers;
    InputFile >> SupportNumbers;
    Support Supports[SupportNumbers];
    for (int i = 0; i < SupportNumbers; ++i)
    {
        Supports[i].ID = i + 1;
        InputFile >> Supports[i].j;
        InputFile >> Supports[i].xR;
        InputFile >> Supports[i].yR;
    }

    // Input Loads
    int LoadNumbers;
    InputFile >> LoadNumbers;
    Load Loads[LoadNumbers];
    for (int i = 0; i < LoadNumbers; ++i)
    {
        Loads[i].ID = i + 1;
        InputFile >> Loads[i].j;
        InputFile >> Loads[i].xC;
        InputFile >> Loads[i].yC;
    }

    // Close Input File
    InputFile.close();

    // Calculate Member's Arguments
    for (int i = 0; i < MemberNumbers; ++i)
    {
        Members[i].l = sqrt(pow((Joints[(Members[i].j - 1)].x - Joints[(Members[i].i) - 1].x), 2) + pow((Joints[(Members[i].j) - 1].y - Joints[(Members[i].i) - 1].y), 2));
        Members[i].Cos = ((Joints[(Members[i].j) - 1].x - Joints[(Members[i].i) -1].x)) / Members[i].l;
        Members[i].Sin = ((Joints[(Members[i].j) - 1].y - Joints[(Members[i].i) -1].y)) / Members[i].l;
        Members[i].AE = Members[i].A * Members[i].E;
    }

    // Make & Calculates Member's K Matrix
    double MembersKMatrix[MemberNumbers][4][4];

    for (int i = 0; i < MemberNumbers; ++i)
    {
        double s2 = pow(Members[i].Sin, 2);
        double c2 = pow(Members[i].Cos, 2);
        double cs = Members[i].Cos * Members[i].Sin;

        MembersKMatrix[i][0][0] = MembersKMatrix[i][2][2] = c2 * Members[i].AE / Members[i].l;
        MembersKMatrix[i][2][0] = MembersKMatrix[i][0][2] = -1 * c2 * Members[i].AE / Members[i].l;
        MembersKMatrix[i][1][1] = MembersKMatrix[i][3][3] = s2 * Members[i].AE / Members[i].l;
        MembersKMatrix[i][1][3] = MembersKMatrix[i][3][1] = -1 * s2 * Members[i].AE / Members[i].l;
        MembersKMatrix[i][0][1] = MembersKMatrix[i][1][0] = MembersKMatrix[i][2][3] = MembersKMatrix[i][3][2] = cs * Members[i].AE / Members[i].l;
        MembersKMatrix[i][2][1] = MembersKMatrix[i][1][2] = MembersKMatrix[i][0][3] = MembersKMatrix[i][3][0] = cs * -1 * Members[i].AE / Members[i].l;
    }

    // Make & Calculates General K Matrix
    double GeneralKMatrix[JointNumbers * 2][JointNumbers * 2];

    for (int i = 0; i < JointNumbers * 2; ++i)
        for (int j = 0; j < JointNumbers * 2; ++j)
            GeneralKMatrix[i][j] = 0;

    for (int i = 0; i < MemberNumbers; ++i)
    {
        int num1 = ((Members[i].i) - 1) * 2;
        int num2 = (Members[i].i) * 2 - 1;  // = num1 + 1
        int num3 = ((Members[i].j) - 1) * 2;
        int num4 = (Members[i].j) * 2 - 1;  // = num3 + 1
        GeneralKMatrix[num1][num1] += MembersKMatrix[i][0][0];
        GeneralKMatrix[num2][num2] += MembersKMatrix[i][1][1];
        GeneralKMatrix[num3][num3] += MembersKMatrix[i][2][2];
        GeneralKMatrix[num4][num4] += MembersKMatrix[i][3][3];
        GeneralKMatrix[num1][num2] += MembersKMatrix[i][0][1];
        GeneralKMatrix[num2][num1] += MembersKMatrix[i][0][1];
        GeneralKMatrix[num1][num3] += MembersKMatrix[i][0][2];
        GeneralKMatrix[num3][num1] += MembersKMatrix[i][0][2];
        GeneralKMatrix[num1][num4] += MembersKMatrix[i][0][3];
        GeneralKMatrix[num4][num1] += MembersKMatrix[i][0][3];
        GeneralKMatrix[num2][num3] += MembersKMatrix[i][1][2];
        GeneralKMatrix[num3][num2] += MembersKMatrix[i][1][2];
        GeneralKMatrix[num2][num4] += MembersKMatrix[i][1][3];
        GeneralKMatrix[num4][num2] += MembersKMatrix[i][1][3];
        GeneralKMatrix[num3][num4] += MembersKMatrix[i][2][3];
        GeneralKMatrix[num4][num3] += MembersKMatrix[i][2][3];
    }

    double PMatrix[JointNumbers * 2];
    for (int i = 0; i < JointNumbers * 2; ++i)
        PMatrix[i] = 0;

    double DeltaMatrix[JointNumbers * 2];
    for (int i = 0; i < JointNumbers * 2; ++i)
        DeltaMatrix[i] = 1368;

    for (int i = 0; i < LoadNumbers; ++i)
    {
        int num1 = (Loads[i].j) * 2 - 2;
        int num2 = (Loads[i].j) * 2 - 1;
        PMatrix[num1] = Loads[i].xC;
        PMatrix[num2] = Loads[i].yC;
    }

    int SNumbers=0;
    for (int i = 0; i < SupportNumbers; ++i)
    {
        int num1 = ((Supports[i].j) - 1) * 2;
        int num2 = (Supports[i].j) * 2 - 1;
        if (Supports[i].xR == 1)
        {
            PMatrix[num1] = 1368;
            DeltaMatrix[num1] = 0;
            SNumbers++;
        }
        if (Supports[i].yR == 1)
        {
            PMatrix[num2]= 1368;
            DeltaMatrix[num2]= 0;
            SNumbers++;
        }
    }

    int FNumbers = JointNumbers * 2 - SNumbers;

    int FNumber[FNumbers];
    int SNumber[SNumbers];

    int counter1 = 0;
    int counter2 = 0;

    for (int i = 0; i < JointNumbers * 2; ++i)
    {
        if (DeltaMatrix[i] == 1368)
        {
            FNumber[counter1] = i;
            counter1++;
        }
        else
        {
            SNumber[counter2] = i;
            counter2++;
        }
    }

    // Make & Calcutate kff,kfs,ksf,kss
    double kff[FNumbers][FNumbers];
    double kss[SNumbers][SNumbers];
    double kfs[FNumbers][SNumbers];
    double ksf[SNumbers][FNumbers];

    for (int i = 0; i < FNumbers; ++i)
        for (int j = 0; j < FNumbers; ++j)
            kff[i][j] = GeneralKMatrix[FNumber[i]][FNumber[j]];
    for (int i = 0; i < FNumbers; ++i)
        for (int j = 0; j < SNumbers; ++j)
            kfs[i][j] = GeneralKMatrix[FNumber[i]][SNumber[j]];
    for (int i = 0; i < SNumbers; ++i)
        for (int j = 0; j < FNumbers; ++j)
            ksf[i][j] = GeneralKMatrix[SNumber[i]][FNumber[j]];
    for (int i = 0; i < SNumbers; ++i)
        for (int j = 0; j < SNumbers; ++j)
            kss[i][j] = GeneralKMatrix[SNumber[i]][SNumber[j]];

    // Make & Calcutate Pf,Ps,Deltaf,Deltas
    double Pf[FNumbers];
    double Ps[SNumbers];
    double Deltaf[FNumbers];
    double Deltas[SNumbers];

    for (int i = 0; i < FNumbers; ++i)
        Pf[i] = PMatrix[FNumber[i]];

    for (int i = 0; i < SNumbers; ++i)
        Ps[i] = PMatrix[SNumber[i]];

    for (int i = 0; i < FNumbers; ++i)
        Deltaf[i] = DeltaMatrix[FNumber[i]];

    for (int i = 0; i < SNumbers; ++i)
        Deltas[i] = DeltaMatrix[SNumber[i]];

    // Make GJMatrix
    double GaussJordanMatrix[FNumbers][FNumbers + 1];
    for (int i = 0; i < FNumbers; ++i)
        for (int j = 0; j < FNumbers; ++j)
            GaussJordanMatrix[i][j] = kff[i][j];

    for (int i = 0; i < FNumbers; ++i)
        GaussJordanMatrix[i][FNumbers]=Pf[i];

    // Solve GJMatrix
    double temp[FNumbers + 1];
    int j = 1;
    while (j < FNumbers)
    {
        for (int i = 0; i < j; ++i)
        {
            for (int k = 0; k < FNumbers + 1; ++k)
                temp[k] = GaussJordanMatrix[j][k] -
                    (GaussJordanMatrix[j][i] / GaussJordanMatrix[i][i]) * GaussJordanMatrix[i][k];
            for (int l = 0; l < FNumbers + 1; ++l)
                GaussJordanMatrix[j][l] = temp[l];
        }
        for (int i = 0; i < j; ++i)
        {
            for (int k = 0; k < FNumbers + 1; ++k)
                temp[k] = GaussJordanMatrix[i][k] -
                    (GaussJordanMatrix[i][j] / GaussJordanMatrix[j][j]) * GaussJordanMatrix[j][k];
            for (int l = 0; l < FNumbers + 1; ++l)
                GaussJordanMatrix[i][l] = temp[l];
        }
        j++;
    }

    // Saving Results of Deltaf
    for (int i = 0; i < FNumbers; ++i)
        Deltaf[i] = GaussJordanMatrix[i][FNumbers] / GaussJordanMatrix[i][i];
    for (int i = 0; i < FNumbers; ++i)
        cout << endl << "Displacement of DOF " << FNumber[i] + 1 << " = " << Deltaf[i] << " cm";
/*
    // Print Delta of Free Joints
    for (int i = 0; i < FNumbers; ++i)
        if (((FNumber[i] + 1) % 2) != 0)
            cout << endl << "  Delta In Joint" << (FNumber[i] + 2) / 2 << ", in X Direction = " << Deltaf[i];
        else
            cout << endl << "  Delta In Joint" << (FNumber[i] + 1) / 2 << ", in Y Direction = " << Deltaf[i];

    // Solve and Make Ps
    double temp2;
    for (int i = 0; i < SNumbers; ++i)
    {
        temp2 = 0;
        for (int j = 0; j < FNumbers; ++j)
            temp2 += ksf[i][j] * Deltaf[j];
        Ps[i] = temp2;
    }

    cout << endl;
    for (int i = 0; i < SNumbers; ++i)
        if (((SNumber[i] + 1) % 2) != 0)
            cout << endl << "  Force In Support" << (SNumber[i] + 2) / 2 << ", in X Direction = " << Ps[i];
        else
            cout << endl << "  Force In Support" << (SNumber[i] + 1) / 2 << ", in Y Direction = " << Ps[i+1];
//*/

    // Making a Matrix for Final Joints
    Joint *FinalJoints;
    FinalJoints = new Joint [JointNumbers];
    for (int i = 0; i < JointNumbers; ++i)
    {
        FinalJoints[i].ID = Joints[i].ID;
        FinalJoints[i].x = Joints[i].x;
        FinalJoints[i].y = Joints[i].y;
    }

    // Makes Final Joints With Deltaf
    for (int i = 0; i < FNumbers; ++i)
        if (((FNumber[i] + 1) % 2) != 0)
            FinalJoints[(FNumber[i] + 2) / 2 - 1].x += Deltaf[i];
        else
            FinalJoints[(FNumber[i] + 1) / 2 - 1].y += Deltaf[i];

    // Make Member's Color
    double e;
    for (int i = 0; i < MemberNumbers; ++i )
    {
        Members[i].l2 = sqrt(pow((FinalJoints[(Members[i].j - 1)].x - FinalJoints[(Members[i].i) - 1].x), 2) + pow((FinalJoints[(Members[i].j) - 1].y - FinalJoints[(Members[i].i) - 1].y), 2));
        Members[i].Stress = (Members[i].l2 - Members[i].l) / Members[i].l * Members[i].E;

        if (Analysed_Truss_Show_Option == "One_Color")
            Members[i].Color = 4; // Aqua
        else if (Analysed_Truss_Show_Option == "Show_Stress")
        {
            e = abs (Members[i].Stress / Members_fy)/5;
            if (e > .6)
                Members[i].Color = 1; // Red
            else if (e > .5 )
                Members[i].Color = 2; // Yellow
            else if (e > .35 )
                Members[i].Color = 73; // Green
            else if (e > 0 )
                Members[i].Color = 3; // Light Green
        }
    }

    MakeDXF(Joints, FinalJoints, Members, Supports, Loads, Original_Truss_Color, 3, 8, 2, JointNumbers, MemberNumbers, SupportNumbers, LoadNumbers);

    getch();

    return 0;
}

void MakeDXF( Joint *Joints, Joint *FinalJoints, Member *Members, Support *Supports,Load *Loads , int MemberColor, int SupportColor, int FinalSupportColor, int LoadColor, int JointCount, int MemberCount, int SupportCount, int LoadCount )
{
    double MaxMemberSize = 0;
    double MaxAE = 0;
    double MaxLoad = 0;

    for (int i = 0; i < MemberCount; ++i)
    {
        if (Members[i].l > MaxMemberSize)
            MaxMemberSize = Members[i].l;
        if ((Members[i].AE) > MaxAE)
            MaxAE = Members[i].AE;
    }
    for (int i = 0; i < LoadCount; ++i)
        Loads[i].Mag = sqrt(Loads[i].xC * Loads[i].xC + Loads[i].yC * Loads[i].yC);
    for (int i = 0; i < LoadCount; ++i)
        if( Loads[i].Mag > MaxLoad )
            MaxLoad = Loads[i].Mag;

    const double SupportSize = MaxMemberSize / 15;
    const double LoadSize = MaxMemberSize / 5 * MaxLoad / MaxAE * 210 * 1.5;
    int iLineCode = 0;
    fstream DXF_Out ;

    DXF_Out.open("Truss Analyser.dxf",ios::out) ;

    DXF_Out << "  2" << endl
            << "HEADER" << endl
            << "  0" << endl
            << "SECTION" << endl
            << "  2" << endl
            << "ENTITIES" << endl;

    // Before Analysis ##################################################################################################################
    for (int i = 0; i < SupportCount; ++i)
        DXF_Out << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << SupportColor << endl
                << " 10" << endl
                << Joints[Supports[i].j - 1].x << endl
                << " 20" << endl
                << Joints[Supports[i].j - 1].y << endl
                << " 11" << endl
                << Joints[Supports[i].j - 1].x - SupportSize << endl
                << " 21" << endl
                << Joints[Supports[i].j - 1].y - SupportSize << endl

                << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << SupportColor << endl
                << " 10" << endl
                << Joints[Supports[i].j - 1].x - SupportSize << endl
                << " 20" << endl
                << Joints[Supports[i].j - 1].y - SupportSize << endl
                << " 11" << endl
                << Joints[Supports[i].j - 1].x + SupportSize << endl
                << " 21" << endl
                << Joints[Supports[i].j - 1].y - SupportSize << endl

                << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << SupportColor << endl
                << " 10" << endl
                << Joints[Supports[i].j - 1].x + SupportSize << endl
                << " 20" << endl
                << Joints[Supports[i].j - 1].y - SupportSize << endl
                << " 11" << endl
                << Joints[Supports[i].j - 1].x << endl
                << " 21" << endl
                << Joints[Supports[i].j - 1].y << endl;

    for ( int i = 0 ; i < MemberCount ; ++i )
        DXF_Out << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << MemberColor << endl
                << " 10" << endl
                << Joints[Members[i].i - 1].x << endl
                << " 20" << endl
                << Joints[Members[i].i - 1].y << endl
                << " 11" << endl
                << Joints[Members[i].j - 1].x << endl
                << " 21" << endl
                << Joints[Members[i].j - 1].y << endl;

    for ( int i = 0 ; i < LoadCount ; ++i )
    {
        if (Loads[i].xC != 0)
        {
        DXF_Out << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << LoadColor << endl
                << " 10" << endl
                << Joints[Loads[i].j - 1].x << endl
                << " 20" << endl
                << Joints[Loads[i].j - 1].y << endl
                << " 11" << endl
                << Joints[Loads[i].j - 1].x - Loads[i].xC / MaxLoad * LoadSize << endl
                << " 21" << endl
                << Joints[Loads[i].j - 1].y << endl

                << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << LoadColor << endl
                << " 10" << endl
                << Joints[Loads[i].j - 1].x << endl
                << " 20" << endl
                << Joints[Loads[i].j - 1].y << endl
                << " 11" << endl
                << Joints[Loads[i].j - 1].x - Loads[i].xC / MaxLoad * LoadSize / 5 << endl
                << " 21" << endl
                << Joints[Loads[i].j - 1].y - Loads[i].xC / MaxLoad * LoadSize / 10 << endl

                << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << LoadColor << endl
                << " 10" << endl
                << Joints[Loads[i].j - 1].x << endl
                << " 20" << endl
                << Joints[Loads[i].j - 1].y << endl
                << " 11" << endl
                << Joints[Loads[i].j - 1].x - Loads[i].xC / MaxLoad * LoadSize/ 5 << endl
                << " 21" << endl
                << Joints[Loads[i].j - 1].y + Loads[i].xC / MaxLoad * LoadSize / 10 << endl;
        }
        if (Loads[i].yC != 0)
        {
        DXF_Out << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << LoadColor << endl
                << " 10" << endl
                << Joints[Loads[i].j - 1].x << endl
                << " 20" << endl
                << Joints[Loads[i].j - 1].y << endl
                << " 11" << endl
                << Joints[Loads[i].j - 1].x << endl
                << " 21" << endl
                << Joints[Loads[i].j - 1].y - Loads[i].yC / MaxLoad * LoadSize << endl

                << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << LoadColor << endl
                << " 10" << endl
                << Joints[Loads[i].j - 1].x << endl
                << " 20" << endl
                << Joints[Loads[i].j - 1].y << endl
                << " 11" << endl
                << Joints[Loads[i].j - 1].x - Loads[i].yC / MaxLoad * LoadSize  / 10 << endl
                << " 21" << endl
                << Joints[Loads[i].j - 1].y - Loads[i].yC / MaxLoad * LoadSize / 5 << endl

                << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << LoadColor << endl
                << " 10" << endl
                << Joints[Loads[i].j - 1].x << endl
                << " 20" << endl
                << Joints[Loads[i].j - 1].y << endl
                << " 11" << endl
                << Joints[Loads[i].j - 1].x + Loads[i].yC / MaxLoad * LoadSize / 10 << endl
                << " 21" << endl
                << Joints[Loads[i].j - 1].y - Loads[i].yC / MaxLoad * LoadSize / 5 << endl;
        }
    }

    // Final Joints ##################################################################################################################
    for ( int i = 0 ; i < SupportCount ; ++i )
        DXF_Out << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << FinalSupportColor << endl
                << " 10" << endl
                << FinalJoints[Supports[i].j - 1].x << endl
                << " 20" << endl
                << FinalJoints[Supports[i].j - 1].y << endl
                << " 11" << endl
                << FinalJoints[Supports[i].j - 1].x - SupportSize << endl
                << " 21" << endl
                << FinalJoints[Supports[i].j - 1].y - SupportSize << endl

                << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << FinalSupportColor << endl
                << " 10" << endl
                << FinalJoints[Supports[i].j - 1].x - SupportSize << endl
                << " 20" << endl
                << FinalJoints[Supports[i].j - 1].y - SupportSize << endl
                << " 11" << endl
                << FinalJoints[Supports[i].j - 1].x + SupportSize << endl
                << " 21" << endl
                << FinalJoints[Supports[i].j - 1].y - SupportSize << endl

                << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << FinalSupportColor << endl
                << " 10" << endl
                << FinalJoints[Supports[i].j - 1].x + SupportSize << endl
                << " 20" << endl
                << FinalJoints[Supports[i].j - 1].y - SupportSize << endl
                << " 11" << endl
                << FinalJoints[Supports[i].j - 1].x << endl
                << " 21" << endl
                << FinalJoints[Supports[i].j - 1].y << endl;

    for ( int i = 0 ; i < MemberCount ; ++i )
        DXF_Out << "  0" << endl
                << "LINE" << endl
                << "  5" << endl
                << ++iLineCode << endl  //Line Code
                << "  8" << endl
                << "0" << endl
                << " 62" << endl
                << "  " << Members[i].Color << endl
                << " 10" << endl
                << FinalJoints[Members[i].i - 1].x << endl
                << " 20" << endl
                << FinalJoints[Members[i].i - 1].y << endl
                << " 11" << endl
                << FinalJoints[Members[i].j - 1].x << endl
                << " 21" << endl
                << FinalJoints[Members[i].j - 1].y << endl;

    DXF_Out << "  0" << endl
            << "ENDSEC" << endl
            << "  0" << endl
            << "EOF" ;
    DXF_Out.close() ;
}
