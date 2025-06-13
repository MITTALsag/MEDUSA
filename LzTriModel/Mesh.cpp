#include "Mesh.h"
#include "MeshTopology.h"
#include "FastDistComp.h"
#include "PatchComputer.h"
#include <LzMath/ToolBox.h>
#include <LzGeom/BBox.h>
#include <LzGeom/Plane3D.h>
#include <LzGeom/Line3D.h>
#include <LzGeom/Point3D.h>
#include <LzGeom/RigidTr3D.h>
#include <LzServices/LzLog.h>
#include <LzServices/Vector.h>
#ifdef USE_CLI
    #include <LzMath/Crypto.h>
    #include <LzServices/XmlReader.h>
    #include <LzServices/XmlWriter.h>
	#include <TxOpenGL/gl.h>
#elif defined(USING_QT) && !defined(NO_OPENGL)
//    #include <LzMath/Crypto.h>
//    #include <LzServices/XmlReader.h>
//    #include <LzServices/XmlWriter.h>
    #include <QtOpenGL> // Need to include this before including glu.h otherwise MSVC does not compile
    #include <GL/glu.h>
#endif
// Qt includes
#ifdef USING_QT
	#include <QtGui/QOpenGLFunctions_2_0>
	#include <QtCore/QTextStream>
	#include <QtCore/QFile>
    #include <QtCore/QFileInfo>
	#include <QDebug>
#endif
#include <fstream>
#include <iostream>
#include <iomanip> // std::precision


//*************** conflicting with std::numeric_limits<double>::max()
#undef max
//*************** conflicting with std::numeric_limits<double>::max()


namespace LzTriModel
{
using LzGeom::Point3D;
using LzGeom::RigidTr3D;
using LzGeom::BBox;
using LzGeom::Line3D;
using LzServices::Vector;
#ifdef USE_CLI
	using System::String;
#endif
using std::ofstream;
//using std::endl; ** conflicts with Qt's *own* endl! Q_CORE_EXPORT QTextStream &endl(QTextStream &s);


#pragma region "Construction & destruction"
//================================================================================
Mesh::Mesh()
{
}

//================================================================================
Mesh::Mesh( const Mesh & pMesh )
{
    *this = pMesh;
}

//================================================================================
Mesh::Mesh( Mesh && pMesh )
{
//    LzLogM("", "*** Mesh::Mesh( Mesh && pMesh )")

    // Geometry
    mVertices = std::move(pMesh.mVertices);
    mTexCoords2D = std::move(pMesh.mTexCoords2D);
    mNormals = std::move(pMesh.mNormals);
    mTriangles = std::move(pMesh.mTriangles);

    // Display lists (a List<uint> is moved here)
    mDLs = std::move(pMesh.mDLs);
}

//================================================================================
Mesh::Mesh( std::function<void(Mesh * pThis)> pInit )
{
    // Call initializer
    pInit(this);
}

//================================================================================
Mesh::~Mesh()
{
    Free();
}

//================================================================================
void Mesh::Free()
{
#if defined(USE_CLI) || defined(USING_QT)
    ResetDisplayList();
#endif

    mVertices.clear();
    mNormals.clear();
    mTriangles.clear();

    mTexCoords2D.Free();
}
#pragma endregion


#pragma region "Load & save"
//================================================================================
void Mesh::Load( const std::string & pFName_Std )
{
	// Read file extension
	std::string lExtension;
#if defined(USING_QT)
	{
		QFileInfo lInfoExt( QString(pFName_Std.c_str()) );
		lExtension = "." + lInfoExt.suffix().toLower().toStdString();
	}
#else
    // Split string into items using "."
    vector<string> lItems;
    LzServices::SplitString( pFName_Std, { "." }, lItems );

    // Check
    if( lItems.size() == 0 )
        LzLogException("", "Cannot extract extension from '"<<pFName_Std<<"'!")

    // Extension = last item
    lExtension = "." + lItems.back();
#endif

    // Dispatch
    if( lExtension == ".xml" )
        LoadXml( pFName_Std );
    else if( lExtension == ".obj" )
        LoadObj( pFName_Std );
    else if( lExtension == ".eobj" )
		LoadEncryptedObj( pFName_Std );
	else if( lExtension == ".mdl" )
        LoadMdl( pFName_Std );
    else if( lExtension == ".stl" )
        LoadStl( pFName_Std );
    else
        LzLogException("", "Unsupported file format: '"<<lExtension<<"'! Cannot load mesh from '"<<pFName_Std<<"'.")

#pragma region "Check consistency"
    try
    {
        for( size_t t=0 ; t<mTriangles.size() ; t++ )
        {
            const Triangle & lTri = mTriangles[t];

            for( int i=0 ; i<3 ; i++ )
            {
                if( lTri.mIdxV[i] >= mVertices.size() )
                    LzLogException("",  "Vertex index out of range (idx= " << lTri.mIdxV[i] << ", max= " << mVertices.size() << ") in triangle " << t << "!" );

                if( mTexCoords2D.Size() )
                {
                    if( lTri.mIdxT[i] >= mTexCoords2D.Size() )
                        LzLogException("",  "Texture vertex index out of range (idx= " << lTri.mIdxT[i] << ", max= " << mTexCoords2D.Size() << ") in triangle " << t << "!" );
                }

                if( lTri.mIdxN[i] >= mNormals.size() )
                    LzLogException("",  "Normal index out of range (idx= " << lTri.mIdxN[i] << ", max= " << mNormals.size() << ") in triangle " << t << "!" );
            }
        }
	}
    catch( const std::exception & e )
	{
        LzLogErr("",  e.what() );
		Free();
        LzLogException("",  "Found inconsistencies in file '" << pFName_Std << "'!" );
	}
#pragma endregion
}

//================================================================================
void Mesh::Save( const std::string & pFName_Std ) const
{
	// Read file extension
	std::string lExtension;
#ifdef USE_CLI
	{
		String ^lSysExtension = System::IO::Path::GetExtension( gcnew String(pFName_Std.c_str()) );
		lSysExtension = lSysExtension->ToLower();
        lExtension = LzServices::StdStrFromCLI( lSysExtension );
	}
#elif defined(USING_QT)
	{
		QFileInfo lInfoExt( QString(pFName_Std.c_str()) );
		lExtension = "." + lInfoExt.suffix().toLower().toStdString();
	}
#else
    // Split string into items using "."
    vector<string> lItems;
    LzServices::SplitString( pFName_Std, { "." }, lItems );

    // Check
    if( lItems.size() == 0 )
        LzLogException("", "Cannot extract extension from '"<<pFName_Std<<"'!")

    // Extension = last item
    lExtension = "." + lItems.back();
#endif

    // Dispatch
    if( lExtension == ".xml" )
        SaveXml( pFName_Std );
    else if( lExtension == ".obj" )
        SaveObj( pFName_Std );
    else if( lExtension == ".eobj" )
		SaveEncryptedObj( pFName_Std );
	else if( lExtension == ".mdl" )
        SaveMdl( pFName_Std );
    else if( lExtension == ".stl" )
        SaveStl( pFName_Std );
    else if( lExtension == ".mesh" )
        SaveMesh( pFName_Std );
    else
       LzLogException("", "Unsupported file format: '"<<lExtension<<"'!");
}

//================================================================================
void Mesh::LoadXml( const std::string & /*pFName_Std*/ )
{
#ifdef USE_CLI
	String ^pFName = gcnew String( pFName_Std.c_str() );

   // Free previous mesh
    Free();

    bool lOk = false;

    _LzLogTry

    // Check & load
    LzServices::XmlReader::CheckXmlFileType( pFName, "TriModel" );
    LzServices::ReadArrayFromXml( pFName, "Vertices", mVertices );
    LzServices::ReadArrayFromXml( pFName, "Normals", mNormals );
    LzServices::ReadArrayFromXml( pFName, "Triangles", mTriangles );

    LzLogMsg("",  "Loaded tri model from: " << pFName_Std << " (" << mVertices.size() << " vertices, " << mNormals.size() << " normals, " << mTriangles.size() << " triangles)." );
    lOk = true;

    _TxLogCatchAndThrow( "Could not read mesh from file: " << pFName_Std << "!" )
    finally
    {
        // If an exception occurred, clean up mesh
        if( !lOk )
            Free();
    }

#else
    LzLogException("", "*** TODO Mesh::LoadXml( const std::string & pFName_Std )");
#endif
}

//================================================================================
void Mesh::SaveXml( const std::string & /*pFName_Std*/ ) const
{
#ifdef USE_CLI
	String ^pFName = gcnew String( pFName_Std.c_str() );

    _LzLogTry

    LzServices::XmlWriter::Open( pFName, "TriModel" );
    LzServices::WriteArrayToXml( "Vertices", mVertices );
    LzServices::WriteArrayToXml( "Normals", mNormals );
    LzServices::WriteArrayToXml( "Triangles", mTriangles );
    LzServices::XmlWriter::Close();

    LzLogMsg("",  "Saved tri model to: " << pFName_Std << " (" << mVertices.size() << " vertices, " << mNormals.size() << " normal(s), " << mTriangles.size() << " triangle(s))." );

    _TxLogCatchAndThrow( "Could not save mesh to file: " << pFName_Std << "!" )

#else
    LzLogException("", "*** TODO Mesh::SaveXml( const std::string & pFName Std) const");
#endif
}


//plus de unsigned int du tout ici²

//faire la chasse a tree3D

//fdc



//================================================================================
void Mesh::LoadEncryptedObj( const std::string & /*pFName_Std*/ )
{
#ifdef USE_CLI
	String ^pFName = gcnew String( pFName_Std.c_str() );

    // Free previous mesh
    Free();

    bool lOk = false;

	System::IO::BinaryReader ^lReader = nullptr;

    _LzLogTry

	lReader = gcnew System::IO::BinaryReader( System::IO::File::Open( pFName, System::IO::FileMode::Open ) );

	// Init cryptosaurus
    LzMath::ToolBox::Crypto lCrypto;
	lCrypto.InitReadFile( lReader );

	//---------------------------------------------
	//
	// NO VALUE CHECKING!
	//
	// Error messages may reveal file structure.
	//
	//---------------------------------------------

	mVertices.resize( lCrypto.Read_UnsInt32(lReader) );
    for( size_t v=0 ; v<mVertices.size() ; v++ )
	{
		mVertices[v].X() = lCrypto.Read_Double64(lReader);
		mVertices[v].Y() = lCrypto.Read_Double64(lReader);
		mVertices[v].Z() = lCrypto.Read_Double64(lReader);
	}

	mTexCoords2D.Resize( lCrypto.Read_UnsInt32(lReader) );
    for( size_t t=0 ; t<mTexCoords2D.Size() ; t++ )
	{
		mTexCoords2D[t].mS = lCrypto.Read_Double64(lReader);
		mTexCoords2D[t].mT = lCrypto.Read_Double64(lReader);
	}

    // Normals
	mNormals.resize( lCrypto.Read_UnsInt32(lReader) );
    for( size_t n=0 ; n<mNormals.size() ; n++ )
	{
		mNormals[n].X() = lCrypto.Read_Double64(lReader);
		mNormals[n].Y() = lCrypto.Read_Double64(lReader);
		mNormals[n].Z() = lCrypto.Read_Double64(lReader);
	}

    // Triangles
	mTriangles.resize( lCrypto.Read_UnsInt32(lReader) );
    for( size_t t=0 ; t<mTriangles.size() ; t++ )
	{
		mTriangles[t].mIdxV[0] = lCrypto.Read_UnsInt32(lReader);
		mTriangles[t].mIdxV[1] = lCrypto.Read_UnsInt32(lReader);
		mTriangles[t].mIdxV[2] = lCrypto.Read_UnsInt32(lReader);

		mTriangles[t].mIdxN[0] = lCrypto.Read_UnsInt32(lReader);
		mTriangles[t].mIdxN[1] = lCrypto.Read_UnsInt32(lReader);
		mTriangles[t].mIdxN[2] = lCrypto.Read_UnsInt32(lReader);

        if( mTexCoords2D.Size() )
		{
			mTriangles[t].mIdxT[0] = lCrypto.Read_UnsInt32(lReader);
			mTriangles[t].mIdxT[1] = lCrypto.Read_UnsInt32(lReader);
			mTriangles[t].mIdxT[2] = lCrypto.Read_UnsInt32(lReader);
		}
	}

     // If no normals loaded
    if( !mNormals.size() )
    {
        // Push default normal
        mNormals.push_back( Vector3D( 1, 0, 0 ) );

        // Log
        LzLogMsg("",  "No normals loaded: pushing default normal: " << mNormals[0].ToString() );
    }

    LzLogMsg("",  "Loaded encrypted tri model from: " << pFName_Std << " (" << mVertices.size() << " vertice(s), " << mTexCoords2D.Size() << " texture coordinate(s), " << mNormals.size() << " normal(s), " << mTriangles.size() << " triangle(s))." );
    lOk = true;

   _TxLogCatchAndThrow( "Could not load encrypted obj mesh from: " << pFName_Std << "!" )
    finally
    {
        if( lReader )
            lReader->Close();

        // If an exception occurred, clean up mesh
        if( !lOk )
            Free();
    }

#else
    LzLogException("", "*** TODO Mesh::LoadEncryptedObj( const std::string & pFName_Std )");
#endif
}

//================================================================================
void Mesh::LoadObj( const std::string & pFName_Std )
{
    // Free previous mesh
    Free();

try
{
    // Open file
    std::ifstream lFile;
    lFile.open( pFName_Std );

    // Check
    if( !lFile.is_open() )
        LzLogException("", "Impossible to open file '"<<pFName_Std<<"'!")

    // Parse file
    string lLine;
    size_t iLine = 0;
    while( std::getline(lFile, lLine) )
    {
        iLine++;

        // Trim line
        LzServices::Trim( lLine );

        // Skip comments
        if( lLine[0] == '#' )
            continue;

        // Skip 'comments'g' (...?)
        if( lLine[0] == 'g' )
            continue;

        // Split line
        vector<string> lTokens;
        LzServices::SplitString( lLine, {" ", "\t"}, lTokens );

        // Skip empty lines
        if( lTokens.size() == 0 )
            continue;

        // Skip material
        if( lTokens[0] == "mtllib"
         || lTokens[0] == "usemtl" )
            continue;

        // Skip groups
        if( lTokens[0] == "g" )
            continue;

        // Read normal
        if( lTokens.size()==4 && lTokens[0]=="vn" )
        {
            // Parse
            Vector3D lNor;
            lNor.X() = LzServices::StrToDouble( lTokens[ 1 ] );
            lNor.Y() = LzServices::StrToDouble( lTokens[ 2 ] );
            lNor.Z() = LzServices::StrToDouble( lTokens[ 3 ] );

            // Stash
            mNormals.push_back( lNor );
        }
        else
        // Read vertex
        if( lTokens.size()==4 && lTokens[0]=="v" )
        {
            // Parse
            Point3D lVer;
            lVer.X() = LzServices::StrToDouble( lTokens[ 1 ] );
            lVer.Y() = LzServices::StrToDouble( lTokens[ 2 ] );
            lVer.Z() = LzServices::StrToDouble( lTokens[ 3 ] );

            // Stash
            mVertices.push_back( lVer );
        }
        else
        // Read texture coordinate
        if( (lTokens.size()==3 || lTokens.size()==4) && lTokens[0]=="vt" )
        {
            // 2D tex coords?
            if( lTokens.size() == 3 )
            {
                // Parse
                TexCoord2D lST;
                lST.mS = LzServices::StrToDouble( lTokens[ 1 ] );
                lST.mT = LzServices::StrToDouble( lTokens[ 2 ] );

                // Stash
                mTexCoords2D.PushBack( lST );
            }
            else
            // 3D tex coords?
            if( lTokens.size() == 3 )
                LzLogException("", "Not implemented! Found 3D texture coordinates in line "<<iLine<<": '"<<lLine<<"'!")
        }
        else
        // Read triangle
        if( (lTokens.size()==4 || lTokens.size()==5) && lTokens[0]=="f" )
        {
            // Triangle?
            if( lTokens.size() == 4 )
            {
                Triangle lTri;
                for( int v=0 ; v<3 ; v++ )
                {
                    // Split each vertex into sub tokens
                    vector<string> lSubTokens;
                    LzServices::SplitString(lTokens[v + 1], {"/"}, lSubTokens, /*KeepEmptyParts*/true);

                    if( lSubTokens.size() == 1 )
                    {
                        // Vertex: "f 1 2 3"
                        lTri.mIdxV[v] = LzServices::StrToUnsInt( lSubTokens[ 0 ] );
                        lTri.mIdxN[v] = 1;
                    }
                    else if( lSubTokens.size() == 2 )
                    {
                        // Vertex/texture: "f 2/2 1/1 3/3"
                        lTri.mIdxV[v] = LzServices::StrToUnsInt( lSubTokens[ 0 ] );
                        lTri.mIdxT[v] = LzServices::StrToUnsInt( lSubTokens[ 1 ] );
                        lTri.mIdxN[v] = 1;
                    }
                    else if( lSubTokens.size() == 3 )
                    {
                        // Vertex/[optional: texture]/normal: "f 2/5/2 1/8/1 3/132/3" or "f 2//2 1//1 3//3"
                        lTri.mIdxV[v] = LzServices::StrToUnsInt( lSubTokens[ 0 ] );
                        if( lSubTokens[ 1 ].length() != 0 )
                            lTri.mIdxT[v] = LzServices::StrToUnsInt( lSubTokens[ 1 ] );
                        lTri.mIdxN[v] = LzServices::StrToUnsInt( lSubTokens[ 2 ] );
                    }
                    else
                        LzLogException("", "Syntax error in line " << iLine << ": '" << lLine << "'!")

                    // Indices in obj start from 1: bring back to 0
                    lTri.mIdxV[v]--;
                    // Update tex coord only if used
                    if( lTri.mIdxT[v] != -1 ) lTri.mIdxT[v]--;
                    lTri.mIdxN[v]--;
                }

                mTriangles.push_back( lTri );
            }
            else
            // Quad?
            if( lTokens.size() == 5 )
                LzLogException("", "Not implemented! Found a quad in line " << iLine << ": '" << lLine << "'!")
        }
        else
            LzLogException("", "Syntax error in line " << iLine << ": '" << lLine << "'!")
    }

    // Push default normal if no normals loaded
    if( !mNormals.size() )
        mNormals.push_back( Vector3D( 1, 0, 0 ) );

    // Log
    LzLogMsg("", "Loaded tri model from: "<<pFName_Std<<" ("<<mVertices.size()<<" vertice(s), "
                 <<mTexCoords2D.Size()<<" texture coordinate(s), "<<mNormals.size()<<" normal(s), "<<mTriangles.size()<<" triangle(s)).")
}
catch( const std::exception & pExc )
{
    // Free possibly partially initialized mesh
    Free();

    // Log lower level exception
    LzLogErr("", "*** Exception: "<<pExc.what())

    // Throw higher level exception
    LzLogException("", "Could not load obj mesh from: "+pFName_Std+"!")
}
}

//================================================================================
void Mesh::SaveEncryptedObj( const std::string & /*pFName_Std*/ ) const
{
#ifdef USE_CLI
	String ^pFName = gcnew String( pFName_Std.c_str() );

	System::IO::BinaryWriter ^lWriter = nullptr;

    _LzLogTry

	lWriter = gcnew System::IO::BinaryWriter( System::IO::File::Open( pFName, System::IO::FileMode::Create ) );

	// Init cryptosaurus
    LzMath::ToolBox::Crypto lCrypto;
	lCrypto.InitWriteFile( lWriter );

    // Vertices
    lCrypto.Write_UnsInt32( mVertices.size(), lWriter );
    for( vector<Point3D>::const_iterator iVer = mVertices.begin() ; iVer != mVertices.end() ; iVer++ )
	{
		lCrypto.Write_Double64( iVer->X(), lWriter );
		lCrypto.Write_Double64( iVer->Y(), lWriter );
		lCrypto.Write_Double64( iVer->Z(), lWriter );
	}

    // Texture coordinates
	lCrypto.Write_UnsInt32( mTexCoords2D.Size(), lWriter );
    for( size_t vt = 0 ; vt < mTexCoords2D.Size() ; vt++ )
	{
		lCrypto.Write_Double64( mTexCoords2D[vt].mS, lWriter );
		lCrypto.Write_Double64( mTexCoords2D[vt].mT, lWriter );
	}

    // Normals
	lCrypto.Write_UnsInt32( mNormals.size(), lWriter );
    for( vector<Vector3D>::const_iterator iNor = mNormals.begin() ; iNor != mNormals.end() ; iNor++ )
	{
		lCrypto.Write_Double64( iNor->X(), lWriter );
		lCrypto.Write_Double64( iNor->Y(), lWriter );
		lCrypto.Write_Double64( iNor->Z(), lWriter );
	}

    // Triangles
	lCrypto.Write_UnsInt32( mTriangles.size(), lWriter );
    for( vector<Triangle>::const_iterator iTri = mTriangles.begin() ; iTri != mTriangles.end() ; iTri++ )
    {
		lCrypto.Write_UnsInt32( iTri->mIdxV[0], lWriter );
		lCrypto.Write_UnsInt32( iTri->mIdxV[1], lWriter );
		lCrypto.Write_UnsInt32( iTri->mIdxV[2], lWriter );

		lCrypto.Write_UnsInt32( iTri->mIdxN[0], lWriter );
		lCrypto.Write_UnsInt32( iTri->mIdxN[1], lWriter );
		lCrypto.Write_UnsInt32( iTri->mIdxN[2], lWriter );

        if( mTexCoords2D.Size() )
		{
			lCrypto.Write_UnsInt32( iTri->mIdxT[0], lWriter );
			lCrypto.Write_UnsInt32( iTri->mIdxT[1], lWriter );
			lCrypto.Write_UnsInt32( iTri->mIdxT[2], lWriter );
		}
    }

    LzLogMsg("",  "Saved encrypted tri model to: " << pFName_Std << " (" << mVertices.size() << " vertice(s), " << mTexCoords2D.Size() << " texture coordinate(s), " << mNormals.size() << " normal(s), " << mTriangles.size() << " triangle(s))." );

    _TxLogCatchAndThrow( "Could not save encrypted obj mesh!" )
    finally
    {
        if( lWriter )
            lWriter->Close();
    }
#else
    LzLogException("", "*** TODO Mesh::SaveEncryptedObj( const std::string & pFName_Std ) const");
#endif
}

//================================================================================
void Mesh::SaveObj( const std::string & pFName_Std ) const
{
    // Open file and check
    ofstream lFile( pFName_Std );
    if( !lFile.is_open() )
        LzLogException("", "Impossible to open file '" << pFName_Std << "'!");

    // Set precision
//lFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    lFile << std::setprecision(15);

    // Vertices
    lFile << "# " << mVertices.size() << " vertice(s)" << std::endl;
    for( vector<Point3D>::const_iterator iVer = mVertices.begin() ; iVer != mVertices.end() ; iVer++ )
        lFile <<  "v " << iVer->X() << " " << iVer->Y() << " " << iVer->Z() << std::endl;
    lFile << std::endl;

    // Texture coordinates
    if( mTexCoords2D.Size() )
    {
        lFile << "# " << mTexCoords2D.Size() << " 2D texture coordinate(s)" << std::endl;
        for( size_t vt=0 ; vt<mTexCoords2D.Size() ; vt++ )
            lFile <<  "vt " << mTexCoords2D[vt].mS << " " << mTexCoords2D[vt].mT << std::endl;
        lFile << std::endl;
    }

    // Normals
    lFile << "# " << mNormals.size() << " normal(s)" << std::endl;
    for( vector<Vector3D>::const_iterator iNor = mNormals.begin() ; iNor != mNormals.end() ; iNor++ )
        lFile << "vn " << iNor->X() << " " << iNor->Y() << " " << iNor->Z() << std::endl;
    lFile << std::endl;

    // Triangles
    lFile << "# " << mTriangles.size() << " triangle(s)" << std::endl;
    for( vector<Triangle>::const_iterator iTri = mTriangles.begin() ; iTri != mTriangles.end() ; iTri++ )
    {
        if( mTexCoords2D.Size() )
            lFile << "f " << (iTri->mIdxV[0] + 1) << "/" << (iTri->mIdxT[0] + 1) << "/" << (iTri->mIdxN[0] + 1) << " "
                          << (iTri->mIdxV[1] + 1) << "/" << (iTri->mIdxT[1] + 1) << "/" << (iTri->mIdxN[1] + 1) << " "
                          << (iTri->mIdxV[2] + 1) << "/" << (iTri->mIdxT[2] + 1) << "/" << (iTri->mIdxN[2] + 1) << std::endl;
        else
            lFile << "f " << (iTri->mIdxV[0] + 1) << "//" << (iTri->mIdxN[0] + 1) << " "
                          << (iTri->mIdxV[1] + 1) << "//" << (iTri->mIdxN[1] + 1) << " "
                          << (iTri->mIdxV[2] + 1) << "//" << (iTri->mIdxN[2] + 1) << std::endl;
    }

    // Log
    LzLogMsg("",  "Saved tri model to: " << pFName_Std << " (" << mVertices.size() << " vertice(s), " << mTexCoords2D.Size() << " texture coordinate(s), " << mNormals.size() << " normal(s), " << mTriangles.size() << " triangle(s))." );
}

//================================================================================
void Mesh::LoadMdl( const std::string & pFName_Std )
{
#ifdef USE_CLI
	String ^pFName = gcnew String( pFName_Std.c_str() );

	// Free previous mesh
    Free();

    bool lOk = false;

    System::IO::TextReader ^lReader = nullptr;

    _LzLogTry

    lReader = gcnew System::IO::StreamReader( pFName );

    cli::array<wchar_t> ^lSpaceTab = { ' ', '\t' };

    // Parse file
    int iLine = 0;
    while( String ^lFileLine = lReader->ReadLine() )
    {
        iLine++;

        //--------------------------------------------------------- Read normals
        if( lFileLine->StartsWith( "[Normals" ) )
        {
            // Already read normals?
            if( mNormals.size() )
                LzLogException("",  "Normals already read from file!" );

            // Read count
            size_t lNorCount = LzServices::StrToUnsInt( lReader->ReadLine() );
            iLine++;

            // Resize and init counter
            mNormals.resize( lNorCount );
            size_t iNorIdx = 0;

            String ^lLine;
            while( ( lLine = lReader->ReadLine() ) && iNorIdx < lNorCount )
            {
                iLine++;

                // Split line
                cli::array<String^> ^lTokens = lLine->Split( lSpaceTab, System::StringSplitOptions::RemoveEmptyEntries );

                // Empty or comment line?
                if( lTokens->Length == 0 || lTokens[0]->StartsWith( "//" ) )
                    continue;

                // Check syntax
                if( lTokens->Length != 3 )
                    LzLogException("", "Could not read normal #"<<iNorIdx<<" in line "<<iLine<<": '"<<LzServices::StdStrFromCLI(lLine)<<"'!");

                Vector3D lNor;
                LzMath::ToolBox::StrToDouble( lTokens[0], lNor.X() );
                LzMath::ToolBox::StrToDouble( lTokens[1], lNor.Y() );
                LzMath::ToolBox::StrToDouble( lTokens[2], lNor.Z() );

                mNormals[iNorIdx++] = lNor;
            }

            // Check underflow
            if( iNorIdx != lNorCount )
                LzLogException("",  "Too few initializers in 'Normals' section!" );
        }
        else
        //--------------------------------------------------------- Read vertices
        if( lFileLine->StartsWith( "[Vertices" ) )
        {
            // Already read vertices?
            if( mVertices.size() )
                LzLogException("",  "Vertices already read from file!" );

            // Read count
            size_t lVerCount = LzServices::StrToUnsInt( lReader->ReadLine() );
            iLine++;

            // Resize and init counter
            mVertices.resize( lVerCount );
            size_t iVerIdx = 0;

            String ^lLine;
            while( ( lLine = lReader->ReadLine() ) && iVerIdx < lVerCount )
            {
                iLine++;

                // Split line
                cli::array<String^> ^lTokens = lLine->Split( lSpaceTab, System::StringSplitOptions::RemoveEmptyEntries );

                // Empty or comment line?
                if( lTokens->Length == 0 || lTokens[0]->StartsWith( "//" ) )
                    continue;

                // Check syntax
                if( lTokens->Length != 3 )
                    LzLogException("", "Could not read vertex #"<<iVerIdx<<" in line "<<iLine<<": '"<<LzServices::StdStrFromCLI(lLine)<<"'!" );

                Point3D lVer;
                LzMath::ToolBox::StrToDouble( lTokens[0], lVer.X() );
                LzMath::ToolBox::StrToDouble( lTokens[1], lVer.Y() );
                LzMath::ToolBox::StrToDouble( lTokens[2], lVer.Z() );

                mVertices[iVerIdx++] = lVer;
            }

            // Check underflow
            if( iVerIdx != lVerCount )
                LzLogException("",  "Too few initializers in 'Vertices' section!" );
        }
        else
        //--------------------------------------------------------- Read triangles
        if( lFileLine->StartsWith( "[Triangles" ) )
        {
            // Already read triangles?
            if( mTriangles.size() )
                LzLogException("",  "Triangles already read from file!" );

            // Read count
            size_t lTriCount = LzServices::StrToUnsInt( lReader->ReadLine() );
            iLine++;

            // Resize and init counter
            mTriangles.resize( lTriCount );
            size_t iTriIdx = 0;

            String ^lLine;
            while( ( lLine = lReader->ReadLine() ) && iTriIdx < lTriCount )
            {
                iLine++;

                // Split line
                cli::array<String^> ^lTokens = lLine->Split( lSpaceTab, System::StringSplitOptions::RemoveEmptyEntries );

                // Empty or comment line?
                if( lTokens->Length == 0 || lTokens[0]->StartsWith( "//" ) )
                    continue;

                // Check syntax
                if( lTokens->Length != 6 )
                    LzLogException("", "Could not read triangle #"<<iTriIdx<<" in line "<<iLine<<": '"<<LzServices::StdStrFromCLI(lLine)<<"'!" );

                Triangle lTri;
                for( int v = 0 ; v < 3 ; v++ )
                {
                    LzMath::ToolBox::StrToUInt32( lTokens[v + 0], lTri.mIdxV[v] );
                    LzMath::ToolBox::StrToUInt32( lTokens[v + 3], lTri.mIdxN[v] );
                }

                mTriangles[iTriIdx++] = lTri;
            }

            // Check underflow
            if( iTriIdx != lTriCount )
                LzLogException("",  "Too few initializers in 'Triangles' section!" );
        }
    }

    // Loaded anything at all?
    if( !mVertices.size() )
        LzLogException("",  "No vertices found in file!" );
    //LzLogErr("", "Warning: no vertices in file!");
    if( !mNormals.size() )
        LzLogErr("",  "Warning: no normals in file!" );
    if( !mTriangles.size() )
        LzLogErr("",  "Warning: no triangles in file!" );

    LzLogMsg("",  "Loaded tri model from: " << pFName_Std << " (" << mVertices.size() << " vertices, " << mNormals.size() << " normals, " << mTriangles.size() << " triangles)." );
    lOk = true;

    _TxLogCatchAndThrow( "Could not load mdl mesh!" )
    finally
    {
        if( lReader )
            lReader->Close();

        // If an exception occurred, clean up mesh
        if( !lOk )
            Free();
    }
#elif defined(USING_QT)
	QFile * lpFile = nullptr;

try
{
    // Free previous mesh
    Free();

	// Open file
    lpFile = new QFile( pFName_Std.c_str() );
    if( !lpFile->open( QFile::ReadOnly ) )
        LzLogException("", "Impossible to open file '" << pFName_Std << "'!")

    QTextStream lReader( lpFile );

	// Parse file
    int iLine = 0;
    QString lLine;
    while( lReader.readLineInto( &lLine ) )
    {
        iLine++;

#pragma region "--------------------------------------------------------- Read normals"
        if( lLine.startsWith( "[Normals" ) )
        {
            // Already read normals?
            if( mNormals.size() )
                LzLogException("", "Normals already read from file!")

            // Read count
			lReader.readLineInto( &lLine );
            iLine++;
            size_t lNorCount;
			{
				bool lParseOk;
				lNorCount = lLine.toUInt( &lParseOk );
				if( !lParseOk )
                    LzLogException("", "Could not read normals count! Error in line "<<iLine<<": '"<<lLine.toStdString()<<"'.")
			}

            // Resize and init counter
            mNormals.resize( lNorCount );
            size_t iNorIdx = 0;

			// Read
            while( lReader.readLineInto(&lLine) && iNorIdx<lNorCount )
            {
                iLine++;

                // Split line
				QStringList lTokens = lLine.split( QRegExp("\\s+"), QString::SkipEmptyParts );

                // Empty or comment line?
                if( lTokens.size()==0 || lTokens[0].startsWith("//") )
                    continue;

                // Check syntax
                if( lTokens.size() != 3 )
                    LzLogException("", "Could not read normal #"<<iNorIdx<<" in line "<<iLine<<": '"<<lLine.toStdString()<<"'!");

				// Lambda
				auto ParseToDouble = [iNorIdx,iLine,lLine]( const QString & pStr, double & pTo )
				{
					bool lParseOk;
					pTo = pStr.toDouble( &lParseOk );
					if( !lParseOk )
                        LzLogException("", "Syntax error in normal #"<<iNorIdx<<" in line "<<iLine<<": '"<<lLine.toStdString()<<"'!");
				};

				// Parse
                Vector3D lNor;
				ParseToDouble( lTokens[0], lNor.X() );
				ParseToDouble( lTokens[1], lNor.Y() );
				ParseToDouble( lTokens[2], lNor.Z() );

				// Stash
                mNormals[iNorIdx++] = lNor;
            }

            // Check underflow
            if( iNorIdx != lNorCount )
                LzLogException("",  "Too few initializers in 'Normals' section!" );
        }
#pragma endregion
        else
#pragma region "--------------------------------------------------------- Read vertices"
        if( lLine.startsWith( "[Vertices" ) )
        {
            // Already read vertices?
            if( mVertices.size() )
                LzLogException("",  "Vertices already read from file!" );

            // Read count
			lReader.readLineInto( &lLine );
            iLine++;
            size_t lVerCount;
			{
				bool lParseOk;
				lVerCount = lLine.toUInt( &lParseOk );
				if( !lParseOk )
                    LzLogException("", "Could not read vertices count! Error in line "<<iLine<<": '"<<lLine.toStdString()<<"'.");
			}

            // Resize and init counter
            mVertices.resize( lVerCount );
            size_t iVerIdx = 0;

            while( lReader.readLineInto(&lLine) && iVerIdx<lVerCount )
            {
                iLine++;

                // Split line
				QStringList lTokens = lLine.split( QRegExp("\\s+"), QString::SkipEmptyParts );

                // Empty or comment line?
                if( lTokens.size()==0 || lTokens[0].startsWith("//") )
                    continue;

                // Check syntax
                if( lTokens.size() != 3 )
                    LzLogException("", "Could not read vertex #"<<iVerIdx<<" in line "<<iLine<<": '"<<lLine.toStdString()<<"'!");

 				// Lambda
				auto ParseToDouble = [iVerIdx,iLine,lLine]( const QString & pStr, double & pTo )
				{
					bool lParseOk;
					pTo = pStr.toDouble( &lParseOk );
					if( !lParseOk )
                        LzLogException("", "Syntax error in vertex #"<<iVerIdx<<" in line "<<iLine<<": '"<<lLine.toStdString()<<"'!");
				};

				// Parse
				Point3D lVer;
				ParseToDouble( lTokens[0], lVer.X() );
				ParseToDouble( lTokens[1], lVer.Y() );
				ParseToDouble( lTokens[2], lVer.Z() );

				// Stash
				mVertices[iVerIdx++] = lVer;
            }

            // Check underflow
            if( iVerIdx != lVerCount )
                LzLogException("",  "Too few initializers in 'Vertices' section!" );
		}
#pragma endregion
		else
#pragma region "--------------------------------------------------------- Read triangles"
        if( lLine.startsWith( "[Triangles" ) )
        {
            // Already read triangles?
            if( mTriangles.size() )
                LzLogException("",  "Triangles already read from file!" );

            // Read count
			lReader.readLineInto( &lLine );
            iLine++;
            size_t lTriCount;
			{
				bool lParseOk;
				lTriCount = lLine.toUInt( &lParseOk );
				if( !lParseOk )
                    LzLogException("", "Could not read triangles count! Error in line "<<iLine<<": '"<<lLine.toStdString()<<"'.");
			}

            // Resize and init counter
            mTriangles.resize( lTriCount );
            size_t iTriIdx = 0;

            while( lReader.readLineInto(&lLine) && iTriIdx<lTriCount )
            {
                iLine++;

				// Split line
				QStringList lTokens = lLine.split( QRegExp("\\s+"), QString::SkipEmptyParts );

                // Empty or comment line?
                if( lTokens.size()==0 || lTokens[0].startsWith("//") )
                    continue;

                // Check syntax
                if( lTokens.size() != 6 )
                    LzLogException("", "Could not read triangle #"<<iTriIdx<<" in line "<<iLine<<": '"<<lLine.toStdString()<<"'!");

 				// Lambda
                auto ParseToUInt = [iTriIdx,iLine,lLine]( const QString & pStr, size_t & pTo )
				{
					bool lParseOk;
					pTo = pStr.toUInt( &lParseOk );
					if( !lParseOk )
                        LzLogException("", "Syntax error in triangle #"<<iTriIdx<<" in line "<<iLine<<": '"<<lLine.toStdString()<<"'!");
				};

				// Parse
                Triangle lTri;
                for( int v = 0 ; v < 3 ; v++ )
                {
                    ParseToUInt( lTokens[v + 0], lTri.mIdxV[v] );
                    ParseToUInt( lTokens[v + 3], lTri.mIdxN[v] );
                }

				// Stash
                mTriangles[iTriIdx++] = lTri;
            }

            // Check underflow
            if( iTriIdx != lTriCount )
                LzLogException("",  "Too few initializers in 'Triangles' section!" );
		}
#pragma endregion
	}

    // Loaded anything at all?
    if( !mVertices.size() )
        LzLogException("",  "No vertices found in file!" );
	//
    if( !mNormals.size() )
        LzLogErr("",  "Warning: no normals in file!" );
	//
    if( !mTriangles.size() )
        LzLogErr("",  "Warning: no triangles in file!" );

	// Log
    LzLogMsg("", "Loaded tri model from: "<<pFName_Std<<" ("<<mVertices.size()<<" vertices, "
             <<mNormals.size()<<" normals, "<<mTriangles.size()<<" triangles).")

	// Close file
	lpFile->close();
	delete lpFile;
}
catch(...)
{
	// Reset object
	Free();

	// Close file if necessary
	if( lpFile )
	{
		if( lpFile->isOpen() )
			lpFile->close();

		delete lpFile;
	}

	// Make some noise
    LzLogException("",  "Could not load mdl mesh!" );
}
#else
    LzLogException("", "*** TODO Mesh::LoadMdl( const std::string & pFName_Std )");
#endif
}

//================================================================================
void Mesh::SaveMdl( const std::string & pFName_Std ) const
{
#ifdef USE_CLI
	String ^pFName = gcnew String( pFName_Std.c_str() );

    System::IO::TextWriter ^lWriter = nullptr;

    _LzLogTry

    lWriter = gcnew System::IO::StreamWriter( pFName );

    // Vertices
    lWriter->WriteLine( "[Vertices, ARRAY1<POINT3D>]" );
    lWriter->WriteLine( "" + mVertices.size() );
    for( vector<Point3D>::const_iterator iVer = mVertices.begin() ; iVer != mVertices.end() ; iVer++ )
        lWriter->WriteLine( "" + iVer->X() + " " + iVer->Y() + " " + iVer->Z() );
    lWriter->WriteLine();

    // Normals
    lWriter->WriteLine( "[Normals, ARRAY1<VECTOR3D>]" );
    lWriter->WriteLine( "" + mNormals.size() );
    for( vector<Vector3D>::const_iterator iNor = mNormals.begin() ; iNor != mNormals.end() ; iNor++ )
        lWriter->WriteLine( "" + iNor->X() + " " + iNor->Y() + " " + iNor->Z() );
    lWriter->WriteLine();

    // Triangles
    lWriter->WriteLine( "[Triangles, ARRAY1<STRING>]" );
    lWriter->WriteLine( "" + mTriangles.size() );
    for( vector<Triangle>::const_iterator iTri = mTriangles.begin() ; iTri != mTriangles.end() ; iTri++ )
        lWriter->WriteLine( "" + iTri->mIdxV[0] + " " + iTri->mIdxV[1] + " " + iTri->mIdxV[2] + " " + iTri->mIdxN[0] + " " + iTri->mIdxN[1] + " " + iTri->mIdxN[2] );

    LzLogMsg("",  "Saved tri model to: " << pFName_Std << " (" << mVertices.size() << " vertices, " << mNormals.size() << " normals, " << mTriangles.size() << " triangles)." );

    _TxLogCatchAndThrow( "Could not save mdl mesh!" )
    finally
    {
        if( lWriter )
            lWriter->Close();
    }

#elif defined(USING_QT)
QFile * lpFile = nullptr;

try
{
	// Open file
	lpFile = new QFile( pFName_Std.c_str() );
	if( !lpFile->open( QFile::WriteOnly ) )
        LzLogException("", "Impossible to open file '" << pFName_Std << "'!");

	// Create writer
	QTextStream lWriter( lpFile );

	// Vertices
    lWriter << "[Vertices, ARRAY1<POINT3D>]" << endl;
    lWriter << mVertices.size() << endl;
    for( size_t v=0 ; v<mVertices.size() ; v++ )
		lWriter << mVertices[v].X() << " " << mVertices[v].Y() << " " << mVertices[v].Z() << endl;
	lWriter << endl;

	// Vertices
	lWriter << "[Normals, ARRAY1<VECTOR3D>]" << endl;
	lWriter << mNormals.size() << endl;
    for( size_t n=0 ; n<mNormals.size() ; n++ )
		lWriter << mNormals[n].X() << " " << mNormals[n].Y() << " " << mNormals[n].Z() << endl;
	lWriter << endl;

	// Triangles
	lWriter << "[Triangles, ARRAY1<STRING>]" << endl;
	lWriter << mTriangles.size() << endl;
    for( size_t t=0 ; t<mTriangles.size() ; t++ )
	{
		const Triangle & lT = mTriangles[t];
		lWriter << lT.mIdxV[0] << " " << lT.mIdxV[1] << " " << lT.mIdxV[2] << " ";
		lWriter << lT.mIdxN[0] << " " << lT.mIdxN[1] << " " << lT.mIdxN[2] << endl;
	}
	lWriter << endl;

    // Log
    LzLogMsg("",  "Saved tri model to: " << pFName_Std << " (" << mVertices.size() << " vertices, " << mNormals.size() << " normals, " << mTriangles.size() << " triangles)." );

	// Clean-up
	lWriter.flush();
	lpFile->close();
	delete lpFile;
}
catch(...)
{
	// Close file if necessary
	if( lpFile )
	{
		if( lpFile->isOpen() )
			lpFile->close();

		delete lpFile;
	}

	// Make some noise
    LzLogException("",  "Could not save mdl mesh!" );
}
#else
    LzLogException("", "*** TODO Mesh::SaveMdl( const std::string & pFName_Std ) const");
#endif
}

//================================================================================
void Mesh::LoadStl( const std::string & /*pFName_Std*/ )
{
#ifdef USE_CLI
	String ^pFName = gcnew String( pFName_Std.c_str() );

    // Free previous mesh
    Free();

    bool lOk = false;

    System::IO::TextReader ^lReader = nullptr;

    _LzLogTry

    lReader = gcnew System::IO::StreamReader( pFName );

    cli::array<wchar_t> ^lSpaceTab = { ' ', '\t' };

#pragma region "******************>>>>>>>>>>>>>>>>>>>>>>>> CODE MINA : achtung baby"

    // Parse file
    int VertexIndex = 0;
    int NormalIndex = 0;
    int iLine = 0;
    String ^lLine;

    while( lLine = lReader->ReadLine() )
    {
        cli::array<String^> ^TempToken = lLine->Split( lSpaceTab, System::StringSplitOptions::RemoveEmptyEntries );

        iLine++;

        if( TempToken->Length == 0 )
            LzLogException("",  "Not recognized syntax in line '"<<iLine<<"'!" )
        else
		if( TempToken[0] == "solid" )
            continue;
        else
		if( TempToken[0] == "endsolid" || ( TempToken->Length == 2 && TempToken[0] == "end" && TempToken[1] == "solid" ) )
            break;
        else
		if( TempToken[0] == "facet" )
        {

            int i;
            Triangle lTri;
            for( i = 0; i < 7; i++ ) // As a facet consists of 7 lines
            {

                cli::array<String^> ^lTokens = lLine->Split( lSpaceTab, System::StringSplitOptions::RemoveEmptyEntries );

                if( lTokens->Length == 0 )
                    LzLogException("", "Not recognized syntax in line '"<<iLine<<"'!");


                // Start by reading the normals
                if( i == 0 )
                {

                    if( lTokens->Length != 5 )
                        LzLogException("", "Not recognized syntax in line '"<<iLine<<"'!");

                    Vector3D lNor;
                    LzMath::ToolBox::StrToDouble( lTokens[2], lNor.X() );
                    LzMath::ToolBox::StrToDouble( lTokens[3], lNor.Y() );
                    LzMath::ToolBox::StrToDouble( lTokens[4], lNor.Z() );
                    mNormals.push_back( lNor );
                    NormalIndex++;
                }


                else if( i == 1 )
                {
                    if( lTokens[0] != "outer" && lTokens[0] != "outerloop" )
                        LzLogException("", "Not recognized syntax in line '"<<iLine<<"'!");
                }

                // Read vertex
                else if( i == 2 || i == 3 || i == 4 )
                {

                    if( lTokens->Length != 4 || lTokens[0] != "vertex" )
                        LzLogException("", "Not recognized syntax in line '"<<iLine<<"'!");

                    Point3D lVer;
                    LzMath::ToolBox::StrToDouble( lTokens[1], lVer.X() );
                    LzMath::ToolBox::StrToDouble( lTokens[2], lVer.Y() );
                    LzMath::ToolBox::StrToDouble( lTokens[3], lVer.Z() );

                    VertexIndex++;
                    mVertices.push_back( lVer );

                }
                else if( i == 5 )
                {
                    if( lTokens[0] != "endloop" && lTokens[0] != "end" )
                        LzLogException("", "Not recognized syntax in line '"<<iLine<<"'!");
                }
                else if( i == 6 )
                {
                    if( lTokens[0] != "endfacet" && lTokens[0] != "end" )
                        LzLogException("", "Not recognized syntax in line '"<<iLine<<"'!");
                }

                if( i != 6 )
                {
                    lLine = lReader->ReadLine();
                    iLine++;
                }
            }
            // Start Constructing the triangles
            lTri.mIdxV[0] = VertexIndex - 3;
            lTri.mIdxV[1] = VertexIndex - 2;
            lTri.mIdxV[2] = VertexIndex - 1;
            lTri.mIdxN[0] = NormalIndex - 1;
            lTri.mIdxN[1] = NormalIndex - 1;
            lTri.mIdxN[2] = NormalIndex - 1;
            mTriangles.push_back( lTri );
        }
        else
            LzLogException("",  "Not recognized syntax in line '"<<iLine<<"'!" );
    }

    LzLogMsg("",  "Number of vertices: " << VertexIndex << "." );
    LzLogMsg("",  "Loaded tri model from: " << pFName_Std << " (" << mVertices.size() << " vertices, " << mNormals.size() << " normals, " << mTriangles.size() << " triangles)." );
    #pragma endregion

    lOk = true;

    _TxLogCatchAndThrow( "Could not load stl mesh from: " << pFName_Std << "!" )
    finally
    {
        if( lReader )
            lReader->Close();

        // If an exception occurred, clean up mesh
        if( !lOk )
            Free();
    }

#else
    LzLogException("", "*** TODO Mesh::LoadStl( const std::string & pFName_Std )");
#endif
}

//================================================================================
void Mesh::SaveStl( const std::string & /*pFName_Std*/ ) const
{
#ifdef USE_CLI
	String ^pFName = gcnew String( pFName_Std.c_str() );

    System::IO::TextWriter ^lWriter = nullptr;

    _LzLogTry

    lWriter = gcnew System::IO::StreamWriter( pFName );

#pragma region "******************>>>>>>>>>>>>>>>>>>>>>>>> CODE MINA : achtung baby"
    // Here the Number of Facet is the number of Triangles
    lWriter->WriteLine( "solid ascii" );

    vector<Point3D>::const_iterator iVer = mVertices.begin();
    vector<Vector3D>::const_iterator iNor = mNormals.begin();
    for( vector<Triangle>::const_iterator iTri = mTriangles.begin(); iTri != mTriangles.end(); iTri++ )
    {
        iNor = mNormals.begin() + iTri->mIdxN[0];
        lWriter->WriteLine( "facet normal:  " + iNor->X() + " " + iNor->Y() + " " + iNor->Z() );

        lWriter->WriteLine( "outer Loop   " );
        //------------------------------------
        iVer = mVertices.begin() + iTri->mIdxV[0];
        lWriter->WriteLine( "vertex    " + iVer->X() + " " + iVer->Y() + " " + iVer->Z() );

        iVer = mVertices.begin() + iTri->mIdxV[1];
        lWriter->WriteLine( "vertex    " + iVer->X() + " " + iVer->Y() + " " + iVer->Z() );

        iVer = mVertices.begin() + iTri->mIdxV[2];
        lWriter->WriteLine( "vertex    " + iVer->X() + " " + iVer->Y() + " " + iVer->Z() );

        //------------------------------------
        lWriter->WriteLine( "end Loop   " );
        lWriter->WriteLine( "end facet   " );
    }

    lWriter->WriteLine( "end solid   " );
#pragma endregion

    LzLogMsg("",  "Saved tri model to: " << pFName_Std << " (" << mVertices.size() << " vertices, " << mNormals.size() << " normals, " << mTriangles.size() << " triangles)." );

    _TxLogCatchAndThrow( "Could not save stl mesh!" )
    finally
    {
        if( lWriter )
            lWriter->Close();
    }

#else
    LzLogException("", "*** TODO Mesh::SaveStl( const std::string & pFName_Std ) const");
#endif
}

//================================================================================
void Mesh::SaveMesh( const std::string & /*pFName_Std*/ ) const
{
#ifdef USE_CLI
	String ^pFName = gcnew String( pFName_Std.c_str() );

    System::IO::TextWriter ^lWriter = nullptr;

    _LzLogTry

    lWriter = gcnew System::IO::StreamWriter( pFName );

    // Header
    lWriter->WriteLine( "MeshVersionFormatted 1" );
    lWriter->WriteLine();
    lWriter->WriteLine( "Dimension" );
    lWriter->WriteLine( "3" );

    // Vertices
    lWriter->WriteLine();
    lWriter->WriteLine( "# Set of mesh vertices" );
    lWriter->WriteLine( "Vertices" );
    lWriter->WriteLine( "" + mVertices.size() );
    for( vector<Point3D>::const_iterator iVer = mVertices.begin() ; iVer != mVertices.end() ; iVer++ )
        lWriter->WriteLine( "" + iVer->X() + " " + iVer->Y() + " " + iVer->Z() + " 0" );

    // Triangles
    lWriter->WriteLine();
    lWriter->WriteLine( "# Set of mesh triangles (v1,v2,v3,tag)" );
    lWriter->WriteLine( "Triangles" );
    lWriter->WriteLine( "" + mTriangles.size() );
    for( vector<Triangle>::const_iterator iTri = mTriangles.begin() ; iTri != mTriangles.end() ; iTri++ )
        lWriter->WriteLine( "" + ( iTri->mIdxV[0] + 1 ) + " " + ( iTri->mIdxV[1] + 1 ) + " " + ( iTri->mIdxV[2] + 1 ) + " 0" );

    // Finish him
    lWriter->WriteLine();
    lWriter->WriteLine( "End" );

    LzLogMsg("",  "Saved tri model to: " << pFName_Std << " (" << mVertices.size() << " vertices, " << mTriangles.size() << " triangles)." );

    _TxLogCatchAndThrow( "Could not save 'mesh' mesh!" )
    finally
    {
        if( lWriter )
            lWriter->Close();
    }

#else
    LzLogException("", "*** TODO Mesh::SaveMesh( const std::string & pFName_Std ) const");
#endif
}
#pragma endregion


#pragma region "Operations"
//================================================================================
bool Mesh::IsMeters( double pMaxLengthInMeters
                     /*=3.0 -- no human body larger than 3 meters...
                     and most body parts are larger than 3 mm*/ ) const
{
    // Check
    if( !IsLoaded() )
        LzLogException("", "Cannot determine scale of an unloaded mesh!")

    // Determine scale
    BBox lBBox( mVertices );
    if( lBBox.Size(lBBox.WidestDim()) < pMaxLengthInMeters )
        return true;
    else
        return false;
}

//================================================================================
size_t Mesh::RemoveUnusedVerticesAndNormals()
{
    // Count unused stuff
    size_t lUnusedVers = 0;
    size_t lUnusedNors = 0;

    // Find unused vertices and normals
    vector<bool> lUsedVer, lUsedNor;
    lUsedVer.assign( mVertices.size(), false );
    lUsedNor.assign( mNormals.size(), false );
    for( size_t t=0 ; t<mTriangles.size() ; t++ )
    {
        const Triangle & lTri = mTriangles[t];

        for( int v = 0 ; v < 3 ; v++ )
        {
            // Fatal error if illegal triangle found
            if( lTri.mIdxV[v] >= mVertices.size() || lTri.mIdxN[v] >= mNormals.size() )
                LzLogException("", "Invalid index found in triangle #"<<t<<", vertex #"<<v<<"! Ver. idx= "<<lTri.mIdxV[v]<<" >= "<<mVertices.size()<<", OR Nor. idx= "<<lTri.mIdxN[v]<<" >= "<<mNormals.size()<<"." )

            if( lTri.mIdxN[v] >= mNormals.size() )
                LzLogException("", "Invalid index found in triangle #"<<t<<", normal #"<<v<<"! Nor. idx= "<<lTri.mIdxN[v]<<" >= "<<mNormals.size()<<"." )

            // Flag vertices and normals as used
            lUsedVer[lTri.mIdxV[v]] = true;
            lUsedNor[lTri.mIdxN[v]] = true;
        }
    }

    // Compute vertex mapping
    vector<int> lOld2NewVers( mVertices.size() );
    {
        vector<Point3D> lNewVers;
        for( size_t v=0 ; v<mVertices.size() ; v++ )
        {
            // Used?
            if( lUsedVer[v] )
            {
                // Store index
                lOld2NewVers[v] = (int)lNewVers.size();

                // Store
                lNewVers.push_back( mVertices[v] );
            }
            else
                lUnusedVers++;
        }

        // Commit
        mVertices = lNewVers;
        LzLogMsg("", "Removed "<<lUnusedVers<<" unused vertices.")
    }

    // Compute normal mapping
    vector<int> lOld2NewNors( mNormals.size() );
    {
        vector<Vector3D> lNewNors;
        for( size_t n = 0 ; n < mNormals.size() ; n++ )
        {
            // Used?
            if( lUsedNor[n] )
            {
                // Store index
                lOld2NewNors[n] = (int)lNewNors.size();

                // Store
                lNewNors.push_back( mNormals[n] );
            }
            else
                lUnusedNors++;
        }

        // Commit
        mNormals = lNewNors;
        LzLogMsg("", "Removed "<<lUnusedNors<<" unused normals.")
    }

    // Re-map triangles
    for( size_t t = 0 ; t < mTriangles.size() ; t++ )
    {
        Triangle & lTri = mTriangles[t];

        for( int v = 0 ; v < 3 ; v++ )
        {
            lTri.mIdxV[v] = lOld2NewVers[ lTri.mIdxV[v] ];
            lTri.mIdxN[v] = lOld2NewNors[ lTri.mIdxN[v] ];
        }
    }

    // Return number of unused items
    return lUnusedVers + lUnusedNors;
}

//================================================================================
void Mesh::RemoveAllTexCoords()
{
    // Clean texture coordinates
    mTexCoords2D.Free();

    // Unref texture coordinates
    for( size_t t = 0 ; t < mTriangles.size() ; t++ )
    {
        for( int i = 0 ; i < 3 ; i++ )
//            mTriangles[t].mIdxT[i] = std::numeric_limits<size_t>::max();
            mTriangles[t].mIdxT[i] = 666;
    }

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
bool Mesh::RemoveDuplicateVertices( double pMergeTol/*=1e-6*/ )
{
    // Compute vertex mapping
    vector<int> lOld2NewVers( mVertices.size() );
    {
        // Find widest dimension
        BBox lBBox( mVertices );
        size_t lBestDim = lBBox.WidestDim();

        // Sort according to widest dimension
        Vector<double> lKeys( mVertices.size() );
        Vector<size_t> lIdx( mVertices.size() );
        for( size_t v = 0 ; v < mVertices.size() ; v++ )
        {
            lKeys[v] = mVertices[v].mV[lBestDim];
            lIdx[v] = v;
        }
        LzMath::ToolBox::QuickSort( lKeys.Buffer(), lKeys.Size(), lIdx.Buffer() );

        // New nodes and mapping
        vector<Point3D> lNewVers;

        // Merge
        size_t lDuplicate = 0;
        size_t iRoot = 0;
        while( iRoot < lIdx.Size() )
        {
            // Stack frame for this class starts at top of the heap
            size_t lClassStart = lNewVers.size();

            // Check all nodes in class
            size_t iNext = iRoot;
            while( iNext < lIdx.Size() )
            {
                // Still in same class?
                if( LzMath::ToolBox::IsZero( lKeys[iRoot] - lKeys[iNext], pMergeTol ) )
                {
                    // Coords of next point in class
                    const Point3D & lNextPos = mVertices[ lIdx[iNext] ];

                    // Try to find identical node
                    int lIdxSame = -1;
                    for( size_t v = lClassStart ; v < lNewVers.size() ; v++ )
                    {
                        // Candidate
                        const Point3D & lPos = lNewVers[v];

                        // Identical?
                        if( LzMath::ToolBox::IsZero( lPos.DistanceTo( lNextPos ), pMergeTol ) )
                        {
                            lIdxSame = v;
                            break;
                        }
                    }

                    // Found identical
                    if( lIdxSame >= 0 )
                    {
                        lDuplicate++;

                        // Link
                        lOld2NewVers[ lIdx[iNext] ] = lIdxSame;
                    }
                    else
                    {
                        // Push new node and link
                        lOld2NewVers[ lIdx[iNext] ] = (int)lNewVers.size();
                        lNewVers.push_back( mVertices[ lIdx[iNext] ] );
                    }
                }
                else
                {
                    // Proceed to next class
                    goto FinitoPepito_Vers;
                }

                // Next
                iNext++;
            }

FinitoPepito_Vers:
            // Start a new class
            iRoot = iNext;
        }

        // Commit
        mVertices = lNewVers;
        LzLogMsg("",  "Merged " << lDuplicate << " duplicate vertices within " << pMergeTol << "." );
    }

    // Compute normal mapping
    //
    // *** nope: not touching the normals
    //

    // Remap triangles
    int lCollapsedCount = 0;
    vector<Triangle> lNewTriangles;
    for( vector<Triangle>::iterator iTri = mTriangles.begin() ; iTri != mTriangles.end() ; iTri++ )
    {
        // Create remapped triangle
        Triangle lNew;
        for( int v = 0 ; v < 3 ; v++ )
        {
            lNew.mIdxV[v] = lOld2NewVers[ iTri->mIdxV[v] ];
            lNew.mIdxT[v] = iTri->mIdxT[v]; //********************************* tex coords should be affected also! (but too complicated)
            lNew.mIdxN[v] = iTri->mIdxN[v];
        }

        // Check for collapsation
        bool lCollapsed = false;
        for( int v = 0 ; v < 3 ; v++ )
        {
            if( lNew.mIdxV[v] == lNew.mIdxV[( v + 1 ) % 3] )
            {
                lCollapsedCount++;
                lCollapsed = true;
                break;
            }
        }

        // Store only if not collapsed
        if( !lCollapsed )
            lNewTriangles.push_back( lNew );
    }

    // Commit
    mTriangles = lNewTriangles;

    // Need to recompute DL
    ResetDisplayList();

    // Log
    if( lCollapsedCount )
    {
        LzLogMsg("",  "Triangles re-mapping: " << lCollapsedCount << " collapsed triangle(s) removed. New vertices might have become unused. Clean again!" );
        return false;
    }
    else
        return true;
}

//================================================================================
bool Mesh::RemoveDuplicateVerticesAndNormals( double pMergeTol/*=1e-6*/ )
{
    // Compute vertex mapping
    vector<int> lOld2NewVers( mVertices.size() );
    {
        // Find widest dimension
        BBox lBBox( mVertices );
        size_t lBestDim = lBBox.WidestDim();

        // Sort according to widest dimension
        Vector<double> lKeys( mVertices.size() );
        Vector<size_t> lIdx( mVertices.size() );
        for( size_t v = 0 ; v < mVertices.size() ; v++ )
        {
            lKeys[v] = mVertices[v].mV[lBestDim];
            lIdx[v] = v;
        }
        LzMath::ToolBox::QuickSort( lKeys.Buffer(), lKeys.Size(), lIdx.Buffer() );

        // New nodes and mapping
        vector<Point3D> lNewVers;

        // Merge
        size_t lDuplicate = 0;
        size_t iRoot = 0;
        while( iRoot < lIdx.Size() )
        {
            // Stack frame for this class starts at top of the heap
            size_t lClassStart = lNewVers.size();

            // Check all nodes in class
            size_t iNext = iRoot;
            while( iNext < lIdx.Size() )
            {
                // Still in same class?
                if( LzMath::ToolBox::IsZero( lKeys[iRoot] - lKeys[iNext], pMergeTol ) )
                {
                    // Coords of next point in class
                    const Point3D & lNextPos = mVertices[ lIdx[iNext] ];

                    // Try to find identical node
                    int lIdxSame = -1;
                    for( size_t v = lClassStart ; v < lNewVers.size() ; v++ )
                    {
                        // Candidate
                        const Point3D & lPos = lNewVers[v];

                        // Identical?
                        if( LzMath::ToolBox::IsZero( lPos.DistanceTo( lNextPos ), pMergeTol ) )
                        {
                            lIdxSame = v;
                            break;
                        }
                    }

                    // Found identical
                    if( lIdxSame >= 0 )
                    {
                        lDuplicate++;

                        // Link
                        lOld2NewVers[ lIdx[iNext] ] = lIdxSame;
                    }
                    else
                    {
                        // Push new node and link
                        lOld2NewVers[ lIdx[iNext] ] = (int)lNewVers.size();
                        lNewVers.push_back( mVertices[ lIdx[iNext] ] );
                    }
                }
                else
                {
                    // Proceed to next class
                    goto FinitoPepito_Vers;
                }

                // Next
                iNext++;
            }

FinitoPepito_Vers:
            // Start a new class
            iRoot = iNext;
        }

        // Commit
        mVertices = lNewVers;
        LzLogMsg("",  "Merged " << lDuplicate << " duplicate vertices within " << pMergeTol << "." );
    }

    // Compute normal mapping
    vector<int> lOld2NewNors( mNormals.size() );
    {
        // Find widest dimension
        BBox lBBox;
        for( size_t n = 0 ; n < mNormals.size() ; n++ )
            lBBox.Update( Point3D( 0, 0, 0 ) + mNormals[n] );
        size_t lBestDim = lBBox.WidestDim();

        // Sort according to widest dimension
        Vector<double> lKeys( mNormals.size() );
        Vector<size_t> lIdx( mNormals.size() );
        for( size_t n = 0 ; n < mNormals.size() ; n++ )
        {
            lKeys[n] = mNormals[n].mV[lBestDim];
            lIdx[n] = n;
        }
        LzMath::ToolBox::QuickSort( lKeys.Buffer(), lKeys.Size(), lIdx.Buffer() );

        // New nodes and mapping
        vector<Vector3D> lNewNors;

        // Merge
        size_t lDuplicate = 0;
        size_t iRoot = 0;
        while( iRoot < lIdx.Size() )
        {
            // Stack frame for this class starts at top of the heap
            size_t lClassStart = lNewNors.size();

            // Check all nodes in class
            size_t iNext = iRoot;
            while( iNext < lIdx.Size() )
            {
                // Still in same class?
                if( LzMath::ToolBox::IsZero( lKeys[iRoot] - lKeys[iNext], pMergeTol ) )
                {
                    // Coords of next normal in class
                    const Vector3D & lNextPos = mNormals[ lIdx[iNext] ];

                    // Try to find identical node
                    int lIdxSame = -1;
                    for( size_t n = lClassStart ; n < lNewNors.size() ; n++ )
                    {
                        // Candidate
                        const Vector3D & lPos = lNewNors[n];

                        // Identical?
                        if( LzMath::ToolBox::IsZero( ( lPos - lNextPos ).Norm(), pMergeTol ) )
                        {
                            lIdxSame = n;
                            break;
                        }
                    }

                    // Found identical
                    if( lIdxSame >= 0 )
                    {
                        lDuplicate++;

                        // Link
                        lOld2NewNors[ lIdx[iNext] ] = lIdxSame;
                    }
                    else
                    {
                        // Push new node and link
                        lOld2NewNors[ lIdx[iNext] ] = (int)lNewNors.size();
                        lNewNors.push_back( mNormals[ lIdx[iNext] ] );
                    }
                }
                else
                {
                    // Proceed to next class
                    goto FinitoPepito_Nors;
                }

                // Next
                iNext++;
            }

FinitoPepito_Nors:
            // Start a new class
            iRoot = iNext;
        }

        // Commit
        mNormals = lNewNors;
        LzLogMsg("",  "Merged " << lDuplicate << " duplicate normals within " << pMergeTol << "." );
    }

    // Remap triangles
    int lCollapsedCount = 0;
    vector<Triangle> lNewTriangles;
    for( vector<Triangle>::iterator iTri = mTriangles.begin() ; iTri != mTriangles.end() ; iTri++ )
    {
        // Create remapped triangle
        Triangle lNew;
        for( int v = 0 ; v < 3 ; v++ )
        {
            lNew.mIdxV[v] = lOld2NewVers[ iTri->mIdxV[v] ];
            lNew.mIdxT[v] = iTri->mIdxT[v];
            lNew.mIdxN[v] = lOld2NewNors[ iTri->mIdxN[v] ];
        }

        // Check for collapsation
        bool lCollapsed = false;
        for( int v = 0 ; v < 3 ; v++ )
        {
            if( lNew.mIdxV[v] == lNew.mIdxV[( v + 1 ) % 3] )
            {
                lCollapsedCount++;
                lCollapsed = true;
                break;
            }
        }

        // Store only if not collapsed
        if( !lCollapsed )
            lNewTriangles.push_back( lNew );
    }

    // Commit
    mTriangles = lNewTriangles;

    // Need to recompute DL
    ResetDisplayList();

    // Log
    if( lCollapsedCount )
    {
        LzLogMsg("",  "Triangles re-mapping: " << lCollapsedCount << " collapsed triangle(s) removed. New vertices might have become unused. Clean again!" );
        return false;
    }
    else
        return true;
}

//================================================================================
size_t Mesh::RemoveCollapsedTriangles()
{
    // Count collapsed triangles
    size_t lCollapsedTris = 0;

    // Check all triangles
    for( size_t t=0 ; t<mTriangles.size() ;)
    {
        // Get triangle
        const Triangle & lTri = mTriangles[t];

        // Check for duplicate indices
        for( int v=0 ; v<3 ; v++ )
        {
            // Same indices
            if( lTri.mIdxV[ v ] == lTri.mIdxV[ (v + 1) % 3 ] )
            {
                // Remove triangle
                mTriangles.erase( mTriangles.begin() + t );

                // Update count
                lCollapsedTris++;

                // Next
                goto NextTriangle;
            }
        }

        // Nope, no pb, increment t
        t++;

    NextTriangle:
        ;
    }

    // Finito pepito
    return lCollapsedTris;
}

//================================================================================
size_t Mesh::RemoveOppositeTriangles()
{
    // Inverse mappings
    static const int sInvMaps[3][3] = {  { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 } };
    //
    auto Opposite = []( const Triangle & pTri0, const Triangle & pTri1 ) -> bool
    {
        // Check all mappings
        for( int map=0 ; map<3 ; map++ )
        {
            if( pTri0.mIdxV[0] == pTri1.mIdxV[ sInvMaps[map][0] ]
             && pTri0.mIdxV[1] == pTri1.mIdxV[ sInvMaps[map][1] ]
             && pTri0.mIdxV[2] == pTri1.mIdxV[ sInvMaps[map][2] ] )
            {
                // Found an inverse mapping
                return true;
            }
        }

        // Didn't find a mapping
        return false;
    };

    // Check all pairs of triangles and count opposite
    size_t lNbOpposite = 0;
    //
    for( size_t t1=0 ; t1<mTriangles.size() ;)
    {
        const Triangle & lTri1 = mTriangles[ t1 ];

        for( size_t t2=t1+1 ; t2<mTriangles.size() ; t2++ )
        {
            const Triangle & lTri2 = mTriangles[ t2 ];

            if( Opposite( lTri1, lTri2 ) )
            {
                // Remove both starting by the one at the end
                mTriangles.erase( mTriangles.begin() + t2 );
                mTriangles.erase( mTriangles.begin() + t1 );

                // Update count
                lNbOpposite += 2;

                // Consider next pair of triangles, keeping t1 where it is
                goto ReTestFromCurr_t1;
            }
        }

        // Next t1
        t1++;

ReTestFromCurr_t1:
        ;
    }

    // Finito pepito
    return lNbOpposite;
}

//================================================================================
Mesh & Mesh::operator=( const Mesh & pMesh )
{
    if( &pMesh != this )
    {
        Free();

        mVertices = pMesh.mVertices;
        mTexCoords2D = pMesh.mTexCoords2D;
        mNormals = pMesh.mNormals;
        mTriangles = pMesh.mTriangles;

        // Leave mDL=0 to force-reset display list
    }

    return *this;
}

//================================================================================
Mesh & Mesh::operator=( Mesh && pMesh )
{
    if( &pMesh != this )
    {
        LzLogM("", "--- MOVING")

        // Geometry
        mVertices = std::move(pMesh.mVertices);
        mTexCoords2D = std::move(pMesh.mTexCoords2D);
        mNormals = std::move(pMesh.mNormals);
        mTriangles = std::move(pMesh.mTriangles);

        // Display lists (a List<uint> is moved here)
        mDLs = std::move(pMesh.mDLs);
    }

    return *this;
}

//================================================================================
void Mesh::Flip( const Plane3D & pMirror )
{
    // Flip vertices
    for( vector<Point3D>::iterator iVer = mVertices.begin() ; iVer != mVertices.end() ; iVer++ )
        *iVer = pMirror.Symmetrical( *iVer );

    // Flip normals
    for( vector<Vector3D>::iterator iNor = mNormals.begin() ; iNor != mNormals.end() ; iNor++ )
        *iNor = pMirror.Symmetrical( *iNor );

    // Flip triangles
    for( vector<Triangle>::iterator iTri = mTriangles.begin() ; iTri != mTriangles.end() ; iTri++ )
        *iTri = Triangle( iTri->mIdxV[0], iTri->mIdxV[2], iTri->mIdxV[1], iTri->mIdxN[0], iTri->mIdxN[2], iTri->mIdxN[1] );

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::Scale( double pScaleX, double pScaleY, double pScaleZ )
{
    for( size_t v = 0 ; v < mVertices.size() ; v++ )
    {
        mVertices[v].X() *= pScaleX;
        mVertices[v].Y() *= pScaleY;
        mVertices[v].Z() *= pScaleZ;
    }

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::RigidTransform( const RigidTr3D & pT )
{
    // Transform vertices
    for( vector<Point3D>::iterator iV = mVertices.begin() ; iV != mVertices.end() ; iV++ )
        ( *iV ) *= pT;

    // Transform normals
    for( vector<Vector3D>::iterator iN = mNormals.begin() ; iN != mNormals.end() ; iN++ )
        ( *iN ) *= pT;

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::AffineTransform( const Matrix & pT )
{
    // Transform vertices
    for( vector<Point3D>::iterator iV=mVertices.begin() ; iV!=mVertices.end() ; iV++ )
        (*iV) *= pT;

    // Compute transform for normals: C = ( (pT)_3x3 )^(-t)
    Matrix A( 3, 3 ), B, C;
    pT.SubMatrixTo( 0, 0, A );
    A.M3x3_InverseTo( B );
    B.TransposeTo( C );

    // Transform normals
    for( vector<Vector3D>::iterator iN=mNormals.begin() ; iN!=mNormals.end() ; iN++ )
    {
        (*iN) *= C;
        iN->Normalize();
    }

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::ElasticTransform( const Vector<Vector3D> & pDispPerVer )
{
	// Check
	if( mVertices.size() != pDispPerVer.Size() )
        LzLogException("", "Invalid number of displacement vectors! "<<mVertices.size()<<" != "<<pDispPerVer.Size()<<".");

	// Apply displacement
    for( size_t v=0 ; v<mVertices.size() ; v++ )
		mVertices[v] += pDispPerVer[v];

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::Merge( const Mesh & pOther )
{
    // Remember initial sizes
    const size_t lInitVerCount = mVertices.size();
    const size_t lInitTexCount = mTexCoords2D.Size();
    const size_t lInitNorCount = mNormals.size();

    // Merge vertices
    for( const Point3D & iVer : pOther.mVertices )
        mVertices.push_back( iVer );

    // Merge texture coordinates
    for( size_t t=0 ; t<pOther.mTexCoords2D.Size() ; t++ )
        mTexCoords2D.PushBack( pOther.mTexCoords2D[t] );

    // Merge normals
    for( const Vector3D & iNor : pOther.mNormals )
        mNormals.push_back( iNor );

    // Merge triangles
    for( const Triangle & iTri : pOther.mTriangles )
    {
        Triangle lTri( lInitVerCount+iTri.mIdxV[0], lInitVerCount+iTri.mIdxV[1], lInitVerCount+iTri.mIdxV[2],
                       lInitNorCount+iTri.mIdxN[0], lInitNorCount+iTri.mIdxN[1], lInitNorCount+iTri.mIdxN[2] );

        // Tex coords: add only if merged mesh had tex coordinates
        if( pOther.mTexCoords2D.Size() )
        {
            for( int i=0 ; i<3 ; i++ )
                lTri.mIdxT[i] = lInitTexCount + iTri.mIdxT[i];
        }

        mTriangles.push_back( lTri );
    }

	// Need to recompute DL
	ResetDisplayList();
}

//================================================================================
void Mesh::MergeTrianglesOnly( const Mesh & pOther )
{
    // Check compatibility
    if( mVertices.size() != pOther.mVertices.size() )
        LzLogException("", "Could not merge triangles only! Vertices size mismatch: "<<mVertices.size()<<" != "<<pOther.mVertices.size()<<".")
    //
    if( mNormals.size() != pOther.mNormals.size() )
        LzLogException("", "Could not merge triangles only! Normals size mismatch: "<<mNormals.size()<<" != "<<pOther.mNormals.size()<<".")
    //
    if( mTexCoords2D.Size() != pOther.mTexCoords2D.Size() )
        LzLogException("", "Could not merge triangles only! Texture coordinates size mismatch: "<<mTexCoords2D.Size()<<" != "<<pOther.mTexCoords2D.Size()<<".")

    // Merge triangles
    for( const Triangle & iTri : pOther.mTriangles )
        mTriangles.push_back( iTri );

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::KeepOnlyTriangles( const vector<size_t> & pTriIdx )
{
    // Keep only the specified triangles
    vector<Triangle> lNewTriangles;
    for( size_t i = 0 ; i < pTriIdx.size() ; i++ )
    {
        size_t lIdx = pTriIdx[i];

        // Check
        if( lIdx >= mTriangles.size() )
            LzLogException("",  "Index " << lIdx << " is out of range! Only " << mTriangles.size() << " triangles in mesh." );

        // Stash
        lNewTriangles.push_back( mTriangles[lIdx] );
    }

    // Commit
    mTriangles = lNewTriangles;

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::SwapTrianglesOrientation()
{
	// Swap triangles
    for( size_t t=0 ; t<mTriangles.size() ; t++ )
	{
		const Triangle & lOldTri = mTriangles[t];
		Triangle lNewTri( lOldTri.mIdxV[2], lOldTri.mIdxV[1], lOldTri.mIdxV[0], lOldTri.mIdxN[2], lOldTri.mIdxN[1], lOldTri.mIdxN[0] );

		mTriangles[t] = lNewTri;
	}

	// Need to recompute DL
	ResetDisplayList();
}

//================================================================================
void Mesh::UglyCut( const Plane3D & pCut )
{
    List<Triangle> lNewTris;

    for( vector<Triangle>::const_iterator iTri = mTriangles.begin() ; iTri != mTriangles.end() ; iTri++ )
    {
        const Point3D & lV0 = mVertices[ iTri->mIdxV[0] ];
        const Point3D & lV1 = mVertices[ iTri->mIdxV[1] ];
        const Point3D & lV2 = mVertices[ iTri->mIdxV[2] ];

        bool lV0_in = pCut.SignedDistanceTo( lV0 ) >= 0;
        bool lV1_in = pCut.SignedDistanceTo( lV1 ) >= 0;
        bool lV2_in = pCut.SignedDistanceTo( lV2 ) >= 0;

        // All in
        if( lV0_in && lV1_in && lV2_in )
            lNewTris.AddTail( *iTri );
    }

    // Store tris
    lNewTris.ToVector( mTriangles );

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::UglyCut( const BBox & pCut )
{
    List<Triangle> lNewTris;

    for( vector<Triangle>::const_iterator iTri = mTriangles.begin() ; iTri != mTriangles.end() ; iTri++ )
    {
        const Point3D & lV0 = mVertices[ iTri->mIdxV[0] ];
        const Point3D & lV1 = mVertices[ iTri->mIdxV[1] ];
        const Point3D & lV2 = mVertices[ iTri->mIdxV[2] ];

        bool lV0_in = pCut.PointIsIn( lV0, true );
        bool lV1_in = pCut.PointIsIn( lV1, true );
        bool lV2_in = pCut.PointIsIn( lV2, true );

        // All in
        if( lV0_in && lV1_in && lV2_in )
            lNewTris.AddTail( *iTri );
    }

    // Store tris
    lNewTris.ToVector( mTriangles );

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
#define V_IDX (lNewVers.Count()+mVertices.size())

//================================================================================
void Mesh::Cut( const Plane3D & pCut, double pSnapRange,
                const vector<Point3D> * pScalpelQuad/*=nullptr*/, const Plane3D * pCutBoundary/*=nullptr*/ )
{
    // Snap vertices
    int snapped = 0;
    for( vector<Point3D>::iterator iVer = mVertices.begin() ; iVer != mVertices.end() ; iVer++ )
    {
//********************** tester les sommets qui snappent dans le quad SI scalpel !!
//********************** ne couper que ces triangles la !
//********************** fermer les outlines qui sortent du scalpel quad

        if( std::abs( pCut.SignedDistanceTo( *iVer ) ) <= pSnapRange )
        {
            *iVer = pCut.Projection( *iVer );
            snapped++;
        }
    }

    // Log
    LzLogMsg("",  "Snapped " << snapped << " vertices within " << pSnapRange << " to cut plane." );

    // Vertices and triangles in new mesh
    List<Point3D> lNewVers;
    List<Triangle> lNewTris;

    for( size_t t=0 ; t<mTriangles.size() ; t++ )
    {
        // Get triangle
        const Triangle & lTri = mTriangles[t];

        // Get vertices
        const Point3D & lV0 = mVertices[ lTri.mIdxV[0] ];
        const Point3D & lV1 = mVertices[ lTri.mIdxV[1] ];
        const Point3D & lV2 = mVertices[ lTri.mIdxV[2] ];

        // Get distances to cut
        const bool lV0_in = pCut.SignedDistanceTo( lV0 ) > -pSnapRange / 2;
        const bool lV1_in = pCut.SignedDistanceTo( lV1 ) > -pSnapRange / 2;
        const bool lV2_in = pCut.SignedDistanceTo( lV2 ) > -pSnapRange / 2;

        // All in thanks to cut boundary
        if( pCutBoundary )
        {
            if( pCutBoundary->SignedDistanceTo(lV0) < 0
             || pCutBoundary->SignedDistanceTo(lV1) < 0
             || pCutBoundary->SignedDistanceTo(lV2) < 0 )
            {
                // One or more vertices are out of bounds ==> triangle must be preserved
                lNewTris.AddTail( lTri );
                continue;
            }
        }

        // All out
        if( !lV0_in && !lV1_in && !lV2_in )
        {
            // Keep triangle only if a scalpel quad is specified
            if( pScalpelQuad )
                lNewTris.AddTail( lTri );

            continue;
        }

        // All in
        if( lV0_in && lV1_in && lV2_in )
        {
            lNewTris.AddTail( lTri );
            continue;
        }

#pragma region "Partial cut"
        char lCutCode = lV0_in << 2 | lV1_in << 1 | lV2_in << 0;
        switch( lCutCode )
        {
        case 1:
        {
            Point3D lI0 = pCut.Intersection( Line3D( lV0, lV2 ) );
            Point3D lI1 = pCut.Intersection( Line3D( lV1, lV2 ) );

            if( pScalpelQuad && !LzMath::ToolBox::SegmentCutsConvexQuad( lI0, lI1, *pScalpelQuad ) )
                lNewTris.AddTail( lTri );
            else
            {
                Triangle lNewT( V_IDX, (V_IDX + 1), lTri.mIdxV[2], lTri.mIdxN[2], lTri.mIdxN[2], lTri.mIdxN[2] );
                lNewTris.AddTail( lNewT );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }
            break;
        }
        case 2:
        {
            Point3D lI0 = pCut.Intersection( Line3D( lV0, lV1 ) );
            Point3D lI1 = pCut.Intersection( Line3D( lV1, lV2 ) );

            if( pScalpelQuad && !LzMath::ToolBox::SegmentCutsConvexQuad( lI0, lI1, *pScalpelQuad ) )
                lNewTris.AddTail( lTri );
            else
            {
                Triangle lNewT( V_IDX, lTri.mIdxV[1], (V_IDX + 1), lTri.mIdxN[1], lTri.mIdxN[1], lTri.mIdxN[1] );
                lNewTris.AddTail( lNewT );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }

            break;
        }
        case 3:
        {
            Point3D lI0 = pCut.Intersection( Line3D( lV0, lV1 ) );
            Point3D lI1 = pCut.Intersection( Line3D( lV0, lV2 ) );

            if( pScalpelQuad && !LzMath::ToolBox::SegmentCutsConvexQuad( lI0, lI1, *pScalpelQuad ) )
                lNewTris.AddTail( lTri );
            else
            {
                Triangle lNewT1(  V_IDX, lTri.mIdxV[1], lTri.mIdxV[2], lTri.mIdxN[1], lTri.mIdxN[1], lTri.mIdxN[2] );
                Triangle lNewT2( (V_IDX + 1),    V_IDX, lTri.mIdxV[2], lTri.mIdxN[2], lTri.mIdxN[1], lTri.mIdxN[2] );
                lNewTris.AddTail( lNewT1 );
                lNewTris.AddTail( lNewT2 );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }

            break;
        }
        case 4:
        {
            Point3D lI0 = pCut.Intersection( Line3D( lV0, lV1 ) );
            Point3D lI1 = pCut.Intersection( Line3D( lV0, lV2 ) );

            if( pScalpelQuad && !LzMath::ToolBox::SegmentCutsConvexQuad( lI0, lI1, *pScalpelQuad ) )
                lNewTris.AddTail( lTri );
            else
            {
                Triangle lNewT( lTri.mIdxV[0], V_IDX, (V_IDX + 1), lTri.mIdxN[0], lTri.mIdxN[0], lTri.mIdxN[0] );
                lNewTris.AddTail( lNewT );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }

            break;
        }
        case 5:
        {
            Point3D lI0 = pCut.Intersection( Line3D( lV0, lV1 ) );
            Point3D lI1 = pCut.Intersection( Line3D( lV1, lV2 ) );

            if( pScalpelQuad && !LzMath::ToolBox::SegmentCutsConvexQuad( lI0, lI1, *pScalpelQuad ) )
                lNewTris.AddTail( lTri );
            else
            {
                Triangle lNewT1( lTri.mIdxV[0], V_IDX,    lTri.mIdxV[2], lTri.mIdxN[0], lTri.mIdxN[0], lTri.mIdxN[2] );
                Triangle lNewT2( lTri.mIdxV[2], V_IDX, (V_IDX + 1), lTri.mIdxN[2], lTri.mIdxN[0], lTri.mIdxN[2] );
                lNewTris.AddTail( lNewT1 );
                lNewTris.AddTail( lNewT2 );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }

            break;
        }
        case 6:
        {
            Point3D lI0 = pCut.Intersection( Line3D( lV0, lV2 ) );
            Point3D lI1 = pCut.Intersection( Line3D( lV1, lV2 ) );

            if( pScalpelQuad && !LzMath::ToolBox::SegmentCutsConvexQuad( lI0, lI1, *pScalpelQuad ) )
                lNewTris.AddTail( lTri );
            else
            {
                Triangle lNewT1( lTri.mIdxV[0], (V_IDX + 1),     V_IDX, lTri.mIdxN[0], lTri.mIdxN[1], lTri.mIdxN[0] );
                Triangle lNewT2( lTri.mIdxV[0],  lTri.mIdxV[1], (V_IDX + 1), lTri.mIdxN[0], lTri.mIdxN[1], lTri.mIdxN[1] );
                lNewTris.AddTail( lNewT1 );
                lNewTris.AddTail( lNewT2 );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }

            break;
        }
        case 0:
        case 7:
        default:
            LzLogException("",  "Illegal cut case!" );
            break;
        }
#pragma endregion
    }

    // Store tris
    lNewTris.ToVector( mTriangles );

    // Store vers
    BrowseList( iV, lNewVers )
    mVertices.push_back( lNewVers.GetAt( iV ) );

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
static Point3D InterpolInterPoint( const Point3D & pA, double pDistA, const Point3D & pB, double pDistB )
{
    if( pDistA * pDistB >= 0 )
        LzLogException("",  "Points are on the same side of the cut! (" << pDistA << " and " << pDistB << ")" );

    return pA + pDistA / ( pDistA - pDistB ) * ( pB - pA );
}

//================================================================================
void Mesh::Cut( const Mesh & pCut, double pSnapRange, CutMeshMode pMode/*, LzTriModel::Patching pPatching=LzTriModel::Patching::Split*/ )
{
    // Set the cutter mesh as a single patch mesh even if it comprises multiple patches
    FastDistComp lFDC;
    lFDC.Set( pCut, LzTriModel::Patching::Split );

    // Snap vertices
    int snapped = 0;
    for( size_t v=0 ; v<mVertices.size() ; v++ )
    {
        Point3D lProjV;
        Vector3D lUnused;
        double lDist = lFDC.SignedDistance( LzGeom::Tree3D::NearestMode::Accurate, mVertices[v], lProjV, lUnused );

        if( std::abs( lDist ) <= pSnapRange )
        {
            mVertices[v] = lProjV;
            snapped++;
        }
    }

    // Log
    LzLogMsg("",  "Snapped " << snapped << " vertices within " << pSnapRange << " to cutting mesh." );

    // Vertices and triangles in new mesh
    List<Point3D> lNewVers;
    List<Triangle> lNewTris;

    for( size_t t = 0 ; t < mTriangles.size() ; t++ )
    {
        const Triangle & lTri = mTriangles[t];

        const Point3D & lV0 = mVertices[ lTri.mIdxV[0] ];
        const Point3D & lV1 = mVertices[ lTri.mIdxV[1] ];
        const Point3D & lV2 = mVertices[ lTri.mIdxV[2] ];

//        Point3D lUnused1;
//        Vector3D lUnused2;
        const double lV0_dist = lFDC.AccurateSignedDistance( lV0 ); //.SignedDistance( LzGeom::Tree3D::NearestMode::Accurate, lV0, lUnused1, lUnused2 );
        const double lV1_dist = lFDC.AccurateSignedDistance( lV1 ); //.SignedDistance( LzGeom::Tree3D::NearestMode::Accurate, lV1, lUnused1, lUnused2 );
        const double lV2_dist = lFDC.AccurateSignedDistance( lV2 ); //.SignedDistance( LzGeom::Tree3D::NearestMode::Accurate, lV2, lUnused1, lUnused2 );
        bool lV0_in = lV0_dist > -pSnapRange / 2;
        bool lV1_in = lV1_dist > -pSnapRange / 2;
        bool lV2_in = lV2_dist > -pSnapRange / 2;

		// Apply remove mode
        if( pMode == CutMeshMode::RemoveOutside )
		{
			// No need to set opposite signed distances
			
			// Need to set opposite flags
			lV0_in = !lV0_in;
			lV1_in = !lV1_in;
			lV2_in = !lV2_in;
		}

        // All out
        if( !lV0_in && !lV1_in && !lV2_in )
            continue;

        // All in
        if( lV0_in && lV1_in && lV2_in )
        {
            lNewTris.AddTail( lTri );
            continue;
        }

#pragma region "Partial cut"
        char lCutCode = lV0_in << 2 | lV1_in << 1 | lV2_in << 0;
        switch( lCutCode )
        {
            case 1:
            {
                Point3D lI0 = InterpolInterPoint( lV0, lV0_dist, lV2, lV2_dist ); //pCut.Intersection( Line3D(lV0,lV2) );
                Point3D lI1 = InterpolInterPoint( lV1, lV1_dist, lV2, lV2_dist ); //pCut.Intersection( Line3D(lV1,lV2) );

                Triangle lNewT( V_IDX, (V_IDX + 1), lTri.mIdxV[2], lTri.mIdxN[2], lTri.mIdxN[2], lTri.mIdxN[2] );
                lNewTris.AddTail( lNewT );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }
            break;

            case 2:
            {
                Point3D lI0 = InterpolInterPoint( lV0, lV0_dist, lV1, lV1_dist ); //pCut.Intersection( Line3D(lV0,lV1) );
                Point3D lI1 = InterpolInterPoint( lV1, lV1_dist, lV2, lV2_dist ); //pCut.Intersection( Line3D(lV1,lV2) );

                Triangle lNewT( V_IDX, lTri.mIdxV[1], (V_IDX + 1), lTri.mIdxN[1], lTri.mIdxN[1], lTri.mIdxN[1] );
                lNewTris.AddTail( lNewT );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }
            break;

            case 3:
            {
                Point3D lI0 = InterpolInterPoint( lV0, lV0_dist, lV1, lV1_dist ); //pCut.Intersection( Line3D(lV0,lV1) );
                Point3D lI1 = InterpolInterPoint( lV0, lV0_dist, lV2, lV2_dist ); //pCut.Intersection( Line3D(lV0,lV2) );

                Triangle lNewT1( V_IDX,   lTri.mIdxV[1], lTri.mIdxV[2], lTri.mIdxN[1], lTri.mIdxN[1], lTri.mIdxN[2] );
                Triangle lNewT2( (V_IDX + 1), V_IDX,          lTri.mIdxV[2], lTri.mIdxN[2], lTri.mIdxN[1], lTri.mIdxN[2] );
                lNewTris.AddTail( lNewT1 );
                lNewTris.AddTail( lNewT2 );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }
            break;

            case 4:
            {
                Point3D lI0 = InterpolInterPoint( lV0, lV0_dist, lV1, lV1_dist ); //pCut.Intersection( Line3D(lV0,lV1) );
                Point3D lI1 = InterpolInterPoint( lV0, lV0_dist, lV2, lV2_dist ); //pCut.Intersection( Line3D(lV0,lV2) );

                Triangle lNewT( lTri.mIdxV[0], V_IDX, (V_IDX + 1), lTri.mIdxN[0], lTri.mIdxN[0], lTri.mIdxN[0] );
                lNewTris.AddTail( lNewT );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }
            break;

            case 5:
            {
                Point3D lI0 = InterpolInterPoint( lV0, lV0_dist, lV1, lV1_dist ); //pCut.Intersection( Line3D(lV0,lV1) );
                Point3D lI1 = InterpolInterPoint( lV1, lV1_dist, lV2, lV2_dist ); //pCut.Intersection( Line3D(lV1,lV2) );

                Triangle lNewT1( lTri.mIdxV[0], V_IDX, lTri.mIdxV[2], lTri.mIdxN[0], lTri.mIdxN[0], lTri.mIdxN[2] );
                Triangle lNewT2( lTri.mIdxV[2], V_IDX, (V_IDX + 1),        lTri.mIdxN[2], lTri.mIdxN[0], lTri.mIdxN[2] );
                lNewTris.AddTail( lNewT1 );
                lNewTris.AddTail( lNewT2 );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }
            break;

            case 6:
            {
                Point3D lI0 = InterpolInterPoint( lV0, lV0_dist, lV2, lV2_dist ); //pCut.Intersection( Line3D(lV0,lV2) );
                Point3D lI1 = InterpolInterPoint( lV1, lV1_dist, lV2, lV2_dist ); //pCut.Intersection( Line3D(lV1,lV2) );

                Triangle lNewT1( lTri.mIdxV[0], (V_IDX + 1),        V_IDX,   lTri.mIdxN[0], lTri.mIdxN[1], lTri.mIdxN[0] );
                Triangle lNewT2( lTri.mIdxV[0], lTri.mIdxV[1], V_IDX + 1, lTri.mIdxN[0], lTri.mIdxN[1], lTri.mIdxN[1] );
                lNewTris.AddTail( lNewT1 );
                lNewTris.AddTail( lNewT2 );

                lNewVers.AddTail( lI0 );
                lNewVers.AddTail( lI1 );
            }
            break;

            case 0:
            case 7:
            default:
                LzLogException("",  "Illegal cut case!" );
                break;
        }
#pragma endregion
    }

    // Store tris
    lNewTris.ToVector( mTriangles );

    // Store vers
    BrowseList( iV, lNewVers )
    mVertices.push_back( lNewVers.GetAt( iV ) );

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::MoveVerticesToMeanPos( size_t pSteps,
                                  const MeshTopology * ppTopo/*=nullptr*/,
                                  const vector<bool> * pIsVerFixed/*=nullptr*/ )
{
    // Pointer to topology: provided by client or my own personal one.
    const MeshTopology * lpTopo;

    MeshTopology lMyOwnTopo;
    if( ppTopo )
        lpTopo = ppTopo;
    else
    {
        lMyOwnTopo.Set( *this );
        lpTopo = &lMyOwnTopo;
    }

    // Link to topology
    const vector<TopVer> & lTopVers = lpTopo->TopVers();
    const vector<TopEdge> & lTopEdges = lpTopo->TopEdges();

    // Check
    if( pIsVerFixed && lTopVers.size() != pIsVerFixed->size() )
        LzLogException("",  "Inconsistent vertex count! Mesh has " << lTopVers.size() << " vertices, while vector has " << pIsVerFixed->size() << "." );

    // Create new table of vertices
    vector<Point3D> lNewVertices( lTopVers.size() );

    // Loop all steps
    for( size_t iStep = 0 ; iStep < pSteps ; iStep++ )
    {
        // Compute mean positions
        for( size_t v = 0 ; v < lTopVers.size() ; v++ )
        {
            // Fixed vertex?
            if( pIsVerFixed && ( *pIsVerFixed )[v] )
            {
                // Don't move!
                lNewVertices[v] = mVertices[v];
            }
            else
            {
                // Edges for this vertex
                const List<size_t> & lEdges = lTopVers[v].mE;

                // New mean point and neighbors count
                Point3D lNewVer = mVertices[v];
                size_t lNeighborsCount = 1 + lEdges.Count();

                BrowseList( iE, lEdges )
                {
                    // Get edge connected vertex
                    const Point3D & lOtherVer = mVertices[ lTopEdges[lEdges.GetAt( iE )].OtherVer( v ) ];

                    // Acc
                    lNewVer.X() += lOtherVer.X();
                    lNewVer.Y() += lOtherVer.Y();
                    lNewVer.Z() += lOtherVer.Z();
                }

                // Normalize
                lNewVer /= lNeighborsCount;

                // Stash
                lNewVertices[v] = lNewVer;
            }
        }

        // Commit
        mVertices = lNewVertices;
    }

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::MoveVerticesToMeanPosAndRecoverBBox( size_t pSteps,
                                                const MeshTopology * ppTopo/*=nullptr*/,
                                                const vector<bool> * pIsVerFixed/*=nullptr*/ )
{
	// Initial bbox
	const BBox lBBox_0( mVertices );

	// Relax mesh
	MoveVerticesToMeanPos( pSteps, ppTopo, pIsVerFixed );

	// Final bbox
	const BBox lBBox_1( mVertices );

	// Check
	if( lBBox_1.Size(lBBox_1.ShortestDim()) < 1e-6 )
	{
		lBBox_1.Log("BBox of mesh after relaxation");
        LzLogException("", "One of the bbox sizes is very very small!");
	}

	// Move new bbox center to origin
	RigidTransform( RigidTr3D(0,0,0, Point3D(0,0,0) - lBBox_1.Center()) );

	// Compensate bbox dims
	Scale( lBBox_0.SizeX()/lBBox_1.SizeX(), lBBox_0.SizeY()/lBBox_1.SizeY(), lBBox_0.SizeZ()/lBBox_1.SizeZ() );

	// Move origin to old bbox center
	RigidTransform( RigidTr3D(0,0,0, lBBox_0.Center() - Point3D(0,0,0)) );
}

//================================================================================
void Mesh::MoveVerticesToMeanPos_HC( size_t pSteps,
                                     const MeshTopology * ppTopo/*=nullptr*/,
                                     const vector<bool> * ppIsVerFixed/*=nullptr*/,
                                     const vector<const Plane3D *> * ppSlidePlanes/*=nullptr*/,
                                     double pAlpha/*=0.4*/, double pBeta/*=0.6*/ )
{
    // Pointer to topology: provided by client or my own personal one.
    const MeshTopology * lpTopo;
	//
    MeshTopology lMyOwnTopo;
    if( ppTopo )
        lpTopo = ppTopo;
    else
    {
        lMyOwnTopo.Set( *this );
        lpTopo = &lMyOwnTopo;
    }

    // Link to topology
    const vector<TopVer> & lTopVers = lpTopo->TopVers();
    const vector<TopEdge> & lTopEdges = lpTopo->TopEdges();

    // Check
    if( ppIsVerFixed && lTopVers.size()!=ppIsVerFixed->size() )
        LzLogException("", "Inconsistent vertex count! Mesh has "<<lTopVers.size()<<" vertices, while fixed vers vector has "<<ppIsVerFixed->size()<<"." );

	// Check
    if( ppSlidePlanes && lTopVers.size()!=ppSlidePlanes->size() )
        LzLogException("", "Inconsistent vertex count! Mesh has "<<lTopVers.size()<<" vertices, while sliding planes vector has "<<ppSlidePlanes->size()<<"." );

	//-------------------------
    // Buffers
	//-------------------------

	// Table of vertices with the original position
	vector<Point3D> Oi;

    // Current configuration Qi is stored in the existing vector mVertices
    vector<Point3D> Pi( lTopVers.size() );  // Mean Values
    Vector<Vector3D> Bi( lTopVers.size() ); // Correction vector to go back towards the original position of v

    // Duplicating the original vertex position in Oi to store it
    if( pSteps > 1 )
        Oi = mVertices;

//#define USE_DEPTH
#ifdef USE_DEPTH
    //const size_t DEPTH = 1; // Equivalent to normal HC Laplacian
    //const size_t DEPTH = 2;
    const size_t DEPTH = 5;

//*********** USING EUCLIDEAN DISTANCE and not topological depth would probably be better...
//*********** USING EUCLIDEAN DISTANCE and not topological depth would probably be better...
//*********** USING EUCLIDEAN DISTANCE and not topological depth would probably be better...

	// Check
	if( DEPTH < 1 )
        LzLogException("", "Depth cannot be less than 1! Please do not feed bullshit params to the algorithms.");

	// Topological neighbours for all vertices
    Vector< List<size_t> > lVerNeighs( lTopVers.size() );
	{
		// For recursive processing
		Vector<bool> lAlreadyPushed;
		lAlreadyPushed.Assign( mVertices.size(), false );

		// Find neighbors for all vertices
		// /!\ A vertex is not listed as one of its neighbors
        for( size_t v=0 ; v<lTopVers.size() ; v++ )
        {
			// Initialize stack with depth=0
            List<size_t> lStack;
			lStack.AddTail( v );
			lAlreadyPushed[v] = true;

			// For all depths
            for( size_t dp=0 ; dp<DEPTH ; dp++ )
			{
				// lStack contains vertices at depth dp

				// Fill lNewStack contains vertices at depth dp
                List<size_t> lNewStack;
				while( lStack.Count() )
				{
					// Pop
                    const size_t w = lStack.GetHead();
					lStack.DelHead();

					// Expand
					const List<int> & lEdges = lTopVers[ w ].mE;
					BrowseList( iE, lEdges )
					{
						// Find opposite ver index
                        const size_t wn = lTopEdges[lEdges.GetAt( iE )].OtherVer( w );

						// Check vertex
						if( !lAlreadyPushed[wn] )
						{
							lNewStack.AddTail( wn );
							lAlreadyPushed[wn] = true;
						}
					}
				}

				// List vertices in new stack
				lVerNeighs[v].Append( lNewStack );

				// Swap stacks if needed
				if( dp + 1 < DEPTH )
					lStack = lNewStack;
			}

			// Clean-up done flags
			lAlreadyPushed[ v ] = false;
            const List<size_t> & lNeighs = lVerNeighs[v];
			BrowseList( iN, lNeighs )
				lAlreadyPushed[ lNeighs.GetAt(iN) ] = false;
		}

//**** DEBUG
//{
//	List<size_t> lSorted;
//	BrowseList( i, lVerNeighs[0] )
//		lSorted.AddIntoIncrList( lVerNeighs[0].GetAt(i), false );
//
//	LzLogMsg("", "Found "<<lSorted.Count()<<" neighbor(s) for vertex 0: "<<lSorted.ListToString());
//	return;
//}
//**** DEBUG

	}
#endif

	//-------------------------
    // Loop all steps
	//-------------------------

    for( size_t iStep=0 ; iStep<pSteps ; iStep++ )
    {
#pragma region "Compute mean positions (smoothing towards the neighbors)"
        for( size_t v=0 ; v<lTopVers.size() ; v++ )
        {
            // Fixed vertex?
            if( ppIsVerFixed && (*ppIsVerFixed)[v] )
            {
                // Don't move!
                Pi[v] = mVertices[v];

//                // Reset B vector
//                Bi[v].Reset();   ===> pointless since the default constructor resets the vector HOWEVER MUST BE UNCOMMENTED IF USING EIGEN
            }
            else
            {
                // Has neighbors?
#ifdef USE_DEPTH
                const List<size_t> & lNeighs = lVerNeighs[v];
                const size_t lNeighborsCount = lNeighs.Count();

                if( lNeighborsCount )
                {
                    // New mean point
                    Pi[v] = Point3D( 0, 0, 0 );

                    BrowseList( iN, lNeighs )
                    {
                        // Get edge connected vertex
                        const Point3D & lOtherVer = mVertices[ lNeighs.GetAt( iN ) ];

                        // Accumulate the position
                        Pi[v].X() += lOtherVer.X();
                        Pi[v].Y() += lOtherVer.Y();
                        Pi[v].Z() += lOtherVer.Z();
                    }

					// Normalize
                    Pi[v] /= lNeighborsCount; // Now we have the mean value Pi (=(all its neighbors)/nbNeigbors)

// Project disp?
if( ppSlidePlanes && (*ppSlidePlanes)[v] )
{
	// Project disp
	Vector3D lDisp = Pi[v] - mVertices[v];
	lDisp = (*ppSlidePlanes)[v]->Projection( lDisp );

	// Apply corrected disp
	Pi[v] = mVertices[v] + lDisp;
}
				}
#else
                const List<size_t> & lEdges = lTopVers[v].mE;
                const size_t lNeighborsCount = lEdges.Count();
				//
                if( lNeighborsCount )
                {
                    // New mean point
                    Pi[v] = Point3D( 0, 0, 0 );

                    BrowseList( iE, lEdges )
                    {
                        // Get edge connected vertex
                        const Point3D & lOtherVer = mVertices[ lTopEdges[lEdges.GetAt( iE )].OtherVer( v ) ];

                        // Accumulate the position
                        Pi[v].X() += lOtherVer.X();
                        Pi[v].Y() += lOtherVer.Y();
                        Pi[v].Z() += lOtherVer.Z();
                    }

					// Normalize
                    Pi[v] /= lNeighborsCount; // Now we have the mean value Pi (=(all its neighbors)/nbNeigbors)

// Project disp?
if( ppSlidePlanes && (*ppSlidePlanes)[v] )
{
	// Project disp
	Vector3D lDisp = Pi[v] - mVertices[v];
	lDisp = (*ppSlidePlanes)[v]->Projection( lDisp );

	// Apply corrected disp
	Pi[v] = mVertices[v] + lDisp;
}
                }
#endif
                else
                {
                    // Don't move!
                    Pi[v] = mVertices[v];
                }

                // Computation of the Bi vector to go back to the original position of v for each edge
                if( pSteps == 1 )       // this if condition could be removed since the else always work, but implies more computation and storage...
                {
                    // This substraction is only valid if only one iteration is applied as the initial position Oi is the same as the previous position qi.
                    // Here Oi and Qi are the same as there is only one iteration, no need to differentiate them and to use alpha
                    Bi[v] = Pi[v] - mVertices[v];
                }
                else
                {
                    // if there is more than one iteration, the stored initial position oi is used as in the paper:  Bi = Pi - (alpha*Oi + (1-alpha)*Qi)
                    Bi[v].X() = Pi[v].X() - ( pAlpha * Oi[v].X() + ( 1 - pAlpha ) * mVertices[v].X() );
                    Bi[v].Y() = Pi[v].Y() - ( pAlpha * Oi[v].Y() + ( 1 - pAlpha ) * mVertices[v].Y() );
                    Bi[v].Z() = Pi[v].Z() - ( pAlpha * Oi[v].Z() + ( 1 - pAlpha ) * mVertices[v].Z() );
                }
            }
        }
#pragma endregion


#pragma region "Computing the vector to go back towards the original position"
        for( size_t v=0 ; v<lTopVers.size() ; v++ )
        {
            // Moving vertex?
            if( !ppIsVerFixed || !(*ppIsVerFixed)[v] )
            {
                // Has neighbors?
                const List<size_t> & lEdges = lTopVers[v].mE;
                size_t lNeighborsCount = lEdges.Count();
                if( lNeighborsCount )
                {
                    // New shift vector
                    Vector3D lNewShift = Vector3D( 0, 0, 0 );

                    BrowseList( iE, lEdges )
                    {
                        // Get vertex connected to vertex v by an edge
                        const Vector3D & lOtherShift = Bi[ lTopEdges[lEdges.GetAt( iE )].OtherVer( v ) ];

                        // Accumulate the vectors Bi (of the neighbors of v) to go back towards the original position of v for each edge
                        lNewShift += lOtherShift;
                    }

                    // Now we have the mean vector Bj weighed by Beta and allowing to go back towards the initial position of v
                    // Finally, Pi = Pi - [ Beta*Bi + (1 - Beta)*(influence of all vectors of v's neighbors) ]
					Vector3D lDisp = pBeta * Bi[v] + ( ( 1 - pBeta ) / lNeighborsCount ) * lNewShift;

// Project disp?
if( ppSlidePlanes && (*ppSlidePlanes)[v] )
	lDisp = (*ppSlidePlanes)[v]->Projection( lDisp );

					// Apply disp
                    Pi[v] -= lDisp;
                }
                else
                {
                    // Don't move!
                    // NOP: Pi[v] is already mVertices[v]
                }
            }
        }
#pragma endregion

        // Commit
        mVertices = Pi;
    }

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::SmoothNormalize( const MeshTopology * ppTopo/*=nullptr*/ )
{
    // Pointer to topology: provided by client or my own personal one.
    const MeshTopology * lpTopo;

    MeshTopology lMyOwnTopo;
    if( ppTopo )
    {
        // Check compatibility
        if( !ppTopo->IsCompatible(*this) )
            LzLogException("", "The provided topology is not compatible with this mesh!")

        // Ok, set
        lpTopo = ppTopo;
    }
    else
    {
        lMyOwnTopo.Set( *this );
        lpTopo = &lMyOwnTopo;
    }

    // Link to topology
    const vector<TopVer> & lTopVers = lpTopo->TopVers();
    const vector<TopTri> & lTopTris = lpTopo->TopTris();

    // Create new table of normals
    vector<Vector3D> lNewNormals( lTopVers.size() );

    // Log counters
    size_t lUncomputableNormals = 0;
    size_t lUnusedVertices = 0;

    // Compute new normals
    for( size_t v = 0 ; v < lTopVers.size() ; v++ )
    {
        // New normal and neighbors count
        Vector3D lNewNor( 0, 0, 0 );
        bool lUsedVertex = false;

        const List<size_t> & lTris = lTopVers[v].mT;
        BrowseList( iT, lTris )
        {
            // Triangle
            const Triangle & lTri = mTriangles[ lTris.GetAt( iT ) ];

            // Vertices
            const Point3D & lV0 = mVertices[ lTri.mIdxV[0] ];
            const Point3D & lV1 = mVertices[ lTri.mIdxV[1] ];
            const Point3D & lV2 = mVertices[ lTri.mIdxV[2] ];

            // Triangle normal
            Vector3D lTriNor = ( lV1 - lV0 ) ^ ( lV2 - lV0 );

            //**********************************************************************************
            // Normal added without prior normalization. As a consequence, the surface of the
            // triangle associated to each TriNor is taken into account.
            //**********************************************************************************

            lNewNor += lTriNor;
            lUsedVertex = true;
        }

        // Normalize
        if( lUsedVertex )
        {
            try
            {
                lNewNor.Normalize();
                lNewNormals[v] = lNewNor;
            }
            catch( ... )
            {
                lUncomputableNormals++;

                // Set default normal
                lNewNormals[v] = Vector3D( 1, 0, 0 );
            }
        }
        else
        {
            lUnusedVertices++;

            // Set default normal
            lNewNormals[v] = Vector3D( 1, 0, 0 );
        }
    }

    // Log
    if( lUncomputableNormals )
        LzLogErr("",  "Could not compute smooth normal for " << lUncomputableNormals << " vertice(s)!" );

    // Log
    if( lUnusedVertices )
        LzLogErr("",  "Found " << lUnusedVertices << " unused vertice(s) in mesh!" );

    // New triangles
    vector<Triangle> lNewTriangles( lTopTris.size() );
    for( size_t t = 0 ; t < lTopTris.size() ; t++ )
    {
        const Triangle & lOldTri = mTriangles[t];
        lNewTriangles[t] = Triangle( lOldTri.mIdxV[0], lOldTri.mIdxV[1], lOldTri.mIdxV[2],
                                     lOldTri.mIdxT[0], lOldTri.mIdxT[1], lOldTri.mIdxT[2],
                                     lOldTri.mIdxV[0], lOldTri.mIdxV[1], lOldTri.mIdxV[2] );
    }

    // Commit
    mNormals = lNewNormals;
    mTriangles = lNewTriangles;

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
static int VerIdxInTri( size_t pVerIdx, const Triangle & pTri )
{
    for( int i=0 ; i<3 ; i++ )
    {
        if( pTri.mIdxV[i] == pVerIdx )
            return i;
    }

    // Mesh topology malfunction...
    LzLogException("",  "Could not find vertex " << pVerIdx << " in triangle " << pTri.ToString() << "! Mesh topology contains an inconsistency..." );
}

//================================================================================
static void FindFanTri( const vector<Triangle> & pTriangles,
                        const vector<TopEdge> & pTopEdges,
                        const vector<TopTri> & pTopTris,
                        const Vector<Vector3D> & pTriNors,
                        const size_t pCurrV,
                        const double pMinSharpAngleDeg,
                        size_t pCurrE,
                        size_t pCurrT,
                        List<size_t> & pTriFan )
{
    while( true )
    {
        const TopEdge & lTopE = pTopEdges[ pCurrE ];

        // Edge is not manifold?
        if( lTopE.mT.Count() != 2 )
            return;

        // Goto next triangle
        size_t lNextT = lTopE.OtherTri( pCurrT );

        // Triangle already listed?
        if( pTriangles[lNextT].mIdxN[ VerIdxInTri( pCurrV, pTriangles[lNextT] ) ] != std::numeric_limits<size_t>::max() )
            return;

        // Jump allowed?
        double lAngleDeg = ( 180.0 / LzMath::PI ) * pTriNors[pCurrT].ShortPositiveAngleTo( pTriNors[lNextT] );
        if( lAngleDeg >= pMinSharpAngleDeg )
            return;

        // Move cursor
        pCurrT = lNextT;
        pCurrE = pTopTris[ pCurrT ].EdgeWithVertex( pCurrV, pCurrE );

        // Stash triangle
        pTriFan.AddTail( pCurrT );
    }
}

//================================================================================
void Mesh::SharpNormalize( double pMinSharpAngleDeg, bool pUseAreaPonderation, const MeshTopology * ppTopo/*=nullptr*/ )
{
try
{
    // Pointer to topology: provided by client or my own personal one.
    const MeshTopology * lpTopo;

    MeshTopology lMyOwnTopo;
    if( ppTopo )
        lpTopo = ppTopo;
    else
    {
        lMyOwnTopo.Set( *this );
        lpTopo = &lMyOwnTopo;
    }


#pragma region "Compute normal vectors and areas for each triangle"
    Vector<double> lTriAreas( mTriangles.size() );
    Vector<Vector3D> lTriNors( mTriangles.size() );
    {
        for( size_t t = 0 ; t < mTriangles.size() ; t++ )
        {
            // The triangle
            const Triangle & lTri = mTriangles[t];

            // The normal
            const Point3D & lV0 = mVertices[ lTri.mIdxV[0] ];
            Vector3D lNor = ( mVertices[ lTri.mIdxV[1] ] - lV0 ) ^ ( mVertices[ lTri.mIdxV[2] ] - lV0 );

            // Store actual area for ponderation, or 1
            lTriAreas[t] = pUseAreaPonderation ? lNor.Norm() / 2.0 : 1 ;

            // Normalize vector
            try
            {
                lNor.Normalize();
                lTriNors[t] = lNor;
            }
            catch( ... )
            {
                // Fake vector with no influence
                lTriNors[t] = Vector3D( 1, 0, 0 );
                lTriAreas[t] = 0;
            }
        }
    }
#pragma endregion


#pragma region "Reset normal mapping in all triangles"
    for( size_t t = 0 ; t < mTriangles.size() ; t++ )
    {
        Triangle & lTri = mTriangles[t];

        for( int i = 0 ; i < 3 ; i++ )
            lTri.mIdxN[i] = std::numeric_limits<size_t>::max();
    }
#pragma endregion


    // New normals (list)
    List<Vector3D> lNewNormalsList;

    // Link to topology
    const vector<TopVer> & lTopVers = lpTopo->TopVers();
    const vector<TopTri> & lTopTris = lpTopo->TopTris();
    const vector<TopEdge> & lTopEdges = lpTopo->TopEdges();


#pragma region "For each vertex, for each triangle around this vertex"
    for( size_t v = 0 ; v < lTopVers.size() ; v++ )
    {
        const TopVer & lTV = lTopVers[ v ];

        // Skip floating vertices
        if( lTV.mT.Count() == 0 )
            continue;

        // Process all triangles within the fan around this vertex
        BrowseList( iT, lTV.mT )
        {
            // The triangle
            const size_t lTriIdx = lTV.mT.GetAt( iT );

            // Find which index has this vertex within the triangle
            const int lVerTriIdx = VerIdxInTri( v, mTriangles[lTriIdx] );

            // This triangle already has a normal at this vertex: skip
            if( mTriangles[lTriIdx].mIdxN[ lVerTriIdx ] != std::numeric_limits<size_t>::max() )
                continue;


#pragma region "Find a smooth triangle fan starting from this triangle, and link vertices to the common fan normal"
            // Link vertex in starting triangle
            mTriangles[lTriIdx].mIdxN[ lVerTriIdx ] = lNewNormalsList.Count();

            // Compute left and right fans
            List<size_t> lTriFanLeft;
            List<size_t> lTriFanRight;
            {
                const TopTri & lTopTri = lTopTris[ lTriIdx ];

                #pragma region "----- Left"
                size_t lEdge0 = lTopTri.EdgeWithVertex( v );
                FindFanTri( mTriangles, lTopEdges, lTopTris, lTriNors, v, pMinSharpAngleDeg, lEdge0, lTriIdx, lTriFanLeft );

                // Link vertices in each triangle of the hemi-fan
                BrowseList( iT, lTriFanLeft )
                {
                    Triangle & lFanTri = mTriangles[ lTriFanLeft.GetAt( iT ) ];
                    lFanTri.mIdxN[ VerIdxInTri( v, lFanTri ) ] = lNewNormalsList.Count();
                }
                #pragma endregion


                #pragma region "----- Right"
                // Go right
                int lEdge1 = lTopTri.EdgeWithVertex( v, lEdge0 );
                FindFanTri( mTriangles, lTopEdges, lTopTris, lTriNors, v, pMinSharpAngleDeg, lEdge1, lTriIdx, lTriFanRight );

                // Link vertices in each triangle of the hemi-fan
                BrowseList( iT, lTriFanRight )
                {
                    Triangle & lFanTri = mTriangles[ lTriFanRight.GetAt( iT ) ];
                    lFanTri.mIdxN[ VerIdxInTri( v, lFanTri ) ] = lNewNormalsList.Count();
                }
                #pragma endregion

////********************************** DEBUG
//{
//
//  LzLogMsg("", "DEBUGGGG*********************************************************");
//
//  List<size_t> lAll;
//  lAll.AddIntoIncrList( lTriIdx, true );
//  BrowseList( iT, lTriFanLeft )
//      lAll.AddIntoIncrList( lTriFanLeft.GetAt( iT ), true );
//  BrowseList( iT, lTriFanRight )
//      lAll.AddIntoIncrList( lTriFanRight.GetAt( iT ), true );
//
//  // Check
//  if( lAll.Count() != lTriFanLeft.Count() + lTriFanRight.Count() + 1 )
//      LzLogException("", "Woooooopsie!!!!!!!");
//}
////********************************** DEBUG
            }
#pragma endregion


            #pragma region "Compute normal for this triangle fan"
            // Init
            Vector3D lFanNor = lTriAreas[ lTriIdx ] * lTriNors[ lTriIdx ];

            // Loop left
            BrowseList( iT, lTriFanLeft )
            {
                size_t lIdx = lTriFanLeft.GetAt( iT );
                lFanNor += lTriAreas[ lIdx ] * lTriNors[ lIdx ];
            }

            // Loop right
            BrowseList( iT, lTriFanRight )
            {
                size_t lIdx = lTriFanRight.GetAt( iT );
                lFanNor += lTriAreas[ lIdx ] * lTriNors[ lIdx ];
            }

            // Normalize
            try
            {
                lFanNor.Normalize();
            }
            catch( ... )
            {
                lFanNor  = Vector3D( 1, 0, 0 );
            }

            // Stash new normal
            lNewNormalsList.AddTail( lFanNor );
            #pragma endregion
        }
    }
#pragma endregion


    // Commit new normals
    lNewNormalsList.ToVector( mNormals );

    // Need to recompute DL
    ResetDisplayList();
}
catch ( const std::exception& e )
{
    LzLogErr("", e.what());
    Free();
    LzLogException("", "Could not sharp normalize mesh!");
}
}

//================================================================================
void Mesh::NormalExtrude( double pNormalShift, bool pUseAreaPonderation )
{
	NormalExtrude( pNormalShift, pNormalShift, pNormalShift, pUseAreaPonderation );
}

//================================================================================
void Mesh::NormalExtrude( double pNormalShift_X, double pNormalShift_Y, double pNormalShift_Z, bool pUseAreaPonderation )
{
    // Compute per-vertex normals
    vector<Vector3D> lPerVertexNormals;
    lPerVertexNormals.assign( mVertices.size(), Vector3D(0, 0, 0) );

    for( vector<Triangle>::const_iterator iTri=mTriangles.begin() ; iTri!=mTriangles.end() ; iTri++ )
    {
        double lTriArea = 1;
        if( pUseAreaPonderation )
        {
            const Point3D & lA = mVertices[ iTri->mIdxV[0] ];
            const Point3D & lB = mVertices[ iTri->mIdxV[1] ];
            const Point3D & lC = mVertices[ iTri->mIdxV[2] ];
            lTriArea = 0.5 * ( ( lB - lA ) ^ ( lC - lA ) ).Norm();
        }

        // Scatter weighted normals
        for( int v=0 ; v<3 ; v++ )
            lPerVertexNormals[ iTri->mIdxV[v] ] += lTriArea * mNormals[ iTri->mIdxN[v] ];
    }

    // Shift vertices
    for( vector<Point3D>::iterator iVer=mVertices.begin() ; iVer!=mVertices.end() ; iVer++ )
    {
        size_t i = iVer - mVertices.begin();
        Vector3D lNormal = lPerVertexNormals[ i ];
        try
        {
            lNormal.Normalize();
        }
        catch( ... )
        {
            // No extrusion here (undefined surface orientation or unused vertex)
            lNormal = Vector3D( 0, 0, 0 );
        }

// Adjust shift vector
lNormal.X() *= pNormalShift_X;
lNormal.Y() *= pNormalShift_Y;
lNormal.Z() *= pNormalShift_Z;

        // Displace vertex
//*iVer = *iVer + pNormalShift * lNormal;
*iVer = *iVer + lNormal;
    }

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::NormalExtrude( const vector<double> & pHeightMap )
{
    // Check
    if( pHeightMap.size() != mVertices.size() )
        LzLogException("", "The provided height map is not compatible with this mesh vertices! (map= "<<pHeightMap.size()<<", vertices= "<<mVertices.size()<<")")

    // Check
    if( mVertices.size() != mNormals.size() )
        LzLogException("", "Cannot normal extrude this mesh! No 1-to-1 mapping between vertices and normals (vertices= "<<mVertices.size()<<", normals= "<<mNormals.size()<<")")

    // Extrude
    for( size_t v=0 ; v<mVertices.size() ; v++ )
    {
        // Displace vertex
        mVertices[v] += pHeightMap[v] * mNormals[ v ];
    }

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::NormalExtrude( const vector<double> & pHeightMap, bool pUseAreaPonderation )
{
    // Check
    if( pHeightMap.size() != mVertices.size() )
        LzLogException("", "The provided height map is not compatible with this mesh vertices! (map= "<<pHeightMap.size()<<", vertices= "<<mVertices.size()<<")")

    // Compute per-vertex normals
    vector<Vector3D> lPerVertexNormals;
    lPerVertexNormals.assign( mVertices.size(), Vector3D(0, 0, 0) );

    for( vector<Triangle>::const_iterator iTri=mTriangles.begin() ; iTri!=mTriangles.end() ; iTri++ )
    {
        double lTriArea = 1;
        if( pUseAreaPonderation )
        {
            const Point3D & lA = mVertices[ iTri->mIdxV[0] ];
            const Point3D & lB = mVertices[ iTri->mIdxV[1] ];
            const Point3D & lC = mVertices[ iTri->mIdxV[2] ];
            lTriArea = 0.5 * ( ( lB - lA ) ^ ( lC - lA ) ).Norm();
        }

        // Scatter weighted normals
        for( int v=0 ; v<3 ; v++ )
            lPerVertexNormals[ iTri->mIdxV[v] ] += lTriArea * mNormals[ iTri->mIdxN[v] ];
    }

    // Shift vertices
    for( size_t v=0 ; v<mVertices.size() ; v++ )
    {
        Point3D & lVer = mVertices[v];
        Vector3D lNormal = lPerVertexNormals[ v ];
        try
        {
            lNormal.Normalize();
        }
        catch( ... )
        {
            // No extrusion here (undefined surface orientation or unused vertex)
            lNormal = Vector3D( 0, 0, 0 );
        }

        // Adjust shift vector
        lNormal *= pHeightMap[v];

        // Displace vertex
        lVer += lNormal;
    }

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::Resample_4xTri()
{
    // Current pointers
    size_t lCurrVerIdx = mVertices.size();
    size_t lCurrTexIdx = mTexCoords2D.Size();
    size_t lCurrNorIdx = mNormals.size();

    vector<Triangle> lNewTriangles;

    for( size_t old_t = 0 ; old_t < mTriangles.size() ; old_t++ )
    {
        const Triangle & lOldTri = mTriangles[old_t];

        // Create new vertices (no references as vectors may be reallocated)
        const Point3D lV0 = mVertices[ lOldTri.mIdxV[0] ];
        const Point3D lV1 = mVertices[ lOldTri.mIdxV[1] ];
        const Point3D lV2 = mVertices[ lOldTri.mIdxV[2] ];

        mVertices.push_back( lV0.MidPointTo( lV1 ) ); // 0 - 1
        mVertices.push_back( lV0.MidPointTo( lV2 ) ); // 0 - 2
        mVertices.push_back( lV1.MidPointTo( lV2 ) ); // 1 - 2

        // Create new normals (no references as vectors may be reallocated)
        const Vector3D lN0 = mNormals[ lOldTri.mIdxN[0] ];
        const Vector3D lN1 = mNormals[ lOldTri.mIdxN[1] ];
        const Vector3D lN2 = mNormals[ lOldTri.mIdxN[2] ];

        Vector3D lN01 = lN0.MeanVectorWith( lN1 );
        Vector3D lN02 = lN0.MeanVectorWith( lN2 );
        Vector3D lN12 = lN1.MeanVectorWith( lN2 );

        try
        {
            lN01.Normalize();
        }
        catch( ... )
        {
            lN01 = lN0; /* because, why not... */
        }
        try
        {
            lN02.Normalize();
        }
        catch( ... )
        {
            lN02 = lN2; /* because, why not... */
        }
        try
        {
            lN12.Normalize();
        }
        catch( ... )
        {
            lN12 = lN1; /* because, why not... */
        }

        mNormals.push_back( lN01 ); // 0 - 1
        mNormals.push_back( lN02 ); // 0 - 2
        mNormals.push_back( lN12 ); // 1 - 2

        // Create new texture coords
        if( lCurrTexIdx )
        {
            const TexCoord2D & lT0 = mTexCoords2D[ lOldTri.mIdxT[0] ];
            const TexCoord2D & lT1 = mTexCoords2D[ lOldTri.mIdxT[1] ];
            const TexCoord2D & lT2 = mTexCoords2D[ lOldTri.mIdxT[2] ];

            TexCoord2D lT01( 0.5 * ( lT0.mS + lT1.mS ), 0.5 * ( lT0.mT + lT1.mT ) );
            TexCoord2D lT02( 0.5 * ( lT0.mS + lT2.mS ), 0.5 * ( lT0.mT + lT2.mT ) );
            TexCoord2D lT12( 0.5 * ( lT1.mS + lT2.mS ), 0.5 * ( lT1.mT + lT2.mT ) );

            mTexCoords2D.PushBack( lT01 ); // 0 - 1
            mTexCoords2D.PushBack( lT02 ); // 0 - 2
            mTexCoords2D.PushBack( lT12 ); // 1 - 2
        }

        // Create triangles
        if( lCurrTexIdx )
        {
            lNewTriangles.push_back( Triangle( lOldTri.mIdxV[0], lCurrVerIdx + 0, lCurrVerIdx + 1, lOldTri.mIdxT[0], lCurrTexIdx + 0, lCurrTexIdx + 1, lOldTri.mIdxN[0], lCurrNorIdx + 0, lCurrNorIdx + 1 ) );
            lNewTriangles.push_back( Triangle( lOldTri.mIdxV[1], lCurrVerIdx + 2, lCurrVerIdx + 0, lOldTri.mIdxT[1], lCurrTexIdx + 2, lCurrTexIdx + 0, lOldTri.mIdxN[1], lCurrNorIdx + 2, lCurrNorIdx + 0 ) );
            lNewTriangles.push_back( Triangle( lOldTri.mIdxV[2], lCurrVerIdx + 1, lCurrVerIdx + 2, lOldTri.mIdxT[2], lCurrTexIdx + 1, lCurrTexIdx + 2, lOldTri.mIdxN[2], lCurrNorIdx + 1, lCurrNorIdx + 2 ) );
            lNewTriangles.push_back( Triangle( lCurrVerIdx + 2, lCurrVerIdx + 1, lCurrVerIdx + 0,    lCurrTexIdx + 2, lCurrTexIdx + 1, lCurrTexIdx + 0,    lCurrNorIdx + 2, lCurrNorIdx + 1, lCurrNorIdx + 0 ) );
        }
        else
        {
            lNewTriangles.push_back( Triangle( lOldTri.mIdxV[0], lCurrVerIdx + 0, lCurrVerIdx + 1, lOldTri.mIdxN[0], lCurrNorIdx + 0, lCurrNorIdx + 1 ) );
            lNewTriangles.push_back( Triangle( lOldTri.mIdxV[1], lCurrVerIdx + 2, lCurrVerIdx + 0, lOldTri.mIdxN[1], lCurrNorIdx + 2, lCurrNorIdx + 0 ) );
            lNewTriangles.push_back( Triangle( lOldTri.mIdxV[2], lCurrVerIdx + 1, lCurrVerIdx + 2, lOldTri.mIdxN[2], lCurrNorIdx + 1, lCurrNorIdx + 2 ) );
            lNewTriangles.push_back( Triangle( lCurrVerIdx + 2, lCurrVerIdx + 1, lCurrVerIdx + 0,    lCurrNorIdx + 2, lCurrNorIdx + 1, lCurrNorIdx + 0 ) );
        }

        // Update pointers
        lCurrVerIdx += 3;
        if( lCurrTexIdx ) lCurrTexIdx += 3;
        lCurrNorIdx += 3;

        //LzLogMsg("", "t= "+old_t);
    }

    // Commit
    mTriangles = lNewTriangles;

    // Need to recompute DL
    ResetDisplayList();
}

//================================================================================
void Mesh::FromBBox( const BBox & pBBox )
{
    // Clean previous
    Free();

    // Vertices
    for( int i = 0; i < 8 ; i++ )
        mVertices.push_back( pBBox.HexaPoint( i ) );

    // Normals
    mNormals.push_back( Vector3D( +1, 0, 0 ) );
    mNormals.push_back( Vector3D( -1, 0, 0 ) );
    mNormals.push_back( Vector3D( 0, +1, 0 ) );
    mNormals.push_back( Vector3D( 0, -1, 0 ) );
    mNormals.push_back( Vector3D( 0, 0, +1 ) );
    mNormals.push_back( Vector3D( 0, 0, -1 ) );

    // Triangles
    mTriangles.push_back( Triangle( 1, 2, 6, 0, 0, 0 ) );
    mTriangles.push_back( Triangle( 6, 5, 1, 0, 0, 0 ) );

    mTriangles.push_back( Triangle( 4, 7, 3, 1, 1, 1 ) );
    mTriangles.push_back( Triangle( 0, 4, 3, 1, 1, 1 ) );

    mTriangles.push_back( Triangle( 2, 3, 6, 2, 2, 2 ) );
    mTriangles.push_back( Triangle( 3, 7, 6, 2, 2, 2 ) );

    mTriangles.push_back( Triangle( 1, 5, 4, 3, 3, 3 ) );
    mTriangles.push_back( Triangle( 1, 4, 0, 3, 3, 3 ) );

    mTriangles.push_back( Triangle( 5, 6, 4, 4, 4, 4 ) );
    mTriangles.push_back( Triangle( 6, 7, 4, 4, 4, 4 ) );

    mTriangles.push_back( Triangle( 3, 1, 0, 5, 5, 5 ) );
    mTriangles.push_back( Triangle( 3, 2, 1, 5, 5, 5 ) );
}

//================================================================================
void Mesh::GetLineIntersections_OLD( const Line3D & pLine, List<Point3D> & pInters,
                                     ScoredPointsList * pScoredInters/*=nullptr*/ ) const
{
    // Remove previous
    pInters.DelAll();

    // Matrices used for local coordinates computation
    Matrix M( 2, 2 ), M_inv( 2, 2 ), B( 2 ), XY( 2 );



// ProstFlex bug: ... but now we need to filter out repeating points!



    // My scored point list
    ScoredPointsList lLocalSortedInters;
    ScoredPointsList * lpSortedInters = pScoredInters ? pScoredInters : &lLocalSortedInters ;

    // Clear previous, if any
    if( pScoredInters )
        pScoredInters->DelAll();



    // Scoring tools
    const Vector3D & lLineDir = pLine.Dir();
    const Point3D & lLineRoot = pLine.Root();

    // Check all triangles
    for( size_t t=0 ; t<mTriangles.size() ; t++ )
    {
        // Triangle and vertices
        const Triangle & lTri = mTriangles[t];
		//
        const Point3D & lA = mVertices[ lTri.mIdxV[0] ];
        const Point3D & lB = mVertices[ lTri.mIdxV[1] ];
        const Point3D & lC = mVertices[ lTri.mIdxV[2] ];

        // Compute intersection
        Point3D lI;
		if( pLine.Intersection_NoExc(Plane3D(lA, lB, lC), lI) )
		{
			// Compute local coordinates
			const Vector3D lX = lB - lA;
			const Vector3D lY = lC - lA;
			const Vector3D lU = lI - lA;

			// Assemble matrix
			M.Elt( 0, 0 ) = lX * lX;
			M.Elt( 1, 1 ) = lY * lY;
			M.Elt( 0, 1 ) = M.Elt( 1, 0 ) = lX * lY;

			// Fill right hand side
			B.Elt( 0 ) = lX * lU;
			B.Elt( 1 ) = lY * lU;

			if( M.M2x2_InverseTo_NoExc(M_inv) )
			{
				M_inv.Mult( B, XY );

				// Check local coordinates
                // ProstFlex: this code is missing some line-mesh intersections...
//                if( XY.Elt( 0 )>=0 && XY.Elt( 1 )>=0 && (XY.Elt( 0 ) + XY.Elt( 1 ))<=1 )
//					pInters.AddTail( lI );
                // ProstFlex bug: Need to add a tolerance margin to capture all entries
                if( XY.Elt( 0 )>=-1e-6 && XY.Elt( 1 )>=-1e-6 && (XY.Elt( 0 ) + XY.Elt( 1 ))<=1+1e-6 )
                    lpSortedInters->AddIntoIncrList( LzServices::Sortable<Point3D,double>(lI, (lI - lLineRoot)*lLineDir), true );
            }
		}
    }

    // ProstFlex bug: ... but now we need to filter out repeating points!
    BrowseList( iSI, *lpSortedInters )
        pInters.AddTail( lpSortedInters->GetAt(iSI).mT );

//LzLogMsg("", "Found "+pInters.Count()+" intersection point(s).");
}

//================================================================================
void Mesh::ComputeTriMats( Vector<Matrix> & pTriMats ) const
{
	// Clear previous
	pTriMats.Free();

	// Compute new
    pTriMats.Resize( mTriangles.size() );
    for( size_t t=0 ; t<mTriangles.size() ; t++ )
    {
        // Triangle and vertices
        const Triangle & lTri = mTriangles[t];
		//
        const Point3D & lA = mVertices[ lTri.mIdxV[0] ];
        const Point3D & lB = mVertices[ lTri.mIdxV[1] ];
        const Point3D & lC = mVertices[ lTri.mIdxV[2] ];

		// Compute local coordinates
		const Vector3D lX = lB - lA;
		const Vector3D lY = lC - lA;

		// Assemble matrix
		Matrix M( 2, 2 );
		M.Elt( 0, 0 ) = lX * lX;
		M.Elt( 1, 1 ) = lY * lY;
		M.Elt( 0, 1 ) = M.Elt( 1, 0 ) = lX * lY;

		// Store inverse of matrix
		Matrix M_inv;
		if( M.M2x2_InverseTo_NoExc(M_inv) )
			pTriMats[ t ] = M_inv;
	}
}

//================================================================================
void Mesh::GetLineIntersections( const Line3D & pLine, List<Point3D> & pInters,
                                 const List<size_t> & pTriIdx, const Vector<Matrix> & pTriMats ) const
{
    // Remove previous
    pInters.DelAll();

    // Matrices used for local coordinates computation
//static Matrix B( 2 ),XY( 2 ); **** 'static' initialization of matrices crashes in debug, in VS2015
Matrix B( 2 ), XY( 2 );

// ProstFlex bug: ... but now we need to filter out repeating points!
List< LzServices::Sortable<Point3D, double> > lSortedInters;
const Vector3D & lLineDir = pLine.Dir();
const Point3D & lLineRoot = pLine.Root();

    // Check selected triangles
BrowseList( iT, pTriIdx )
//for( size_t t=0 ; t<mTriangles.size() ; t++ )
   {
		// Read index
const size_t lIdx = pTriIdx.GetAt( iT );
//const size_t lIdx = t;

        // Triangle and vertices
        const Triangle & lTri = mTriangles[ lIdx ];
		//
        const Point3D & lA = mVertices[ lTri.mIdxV[0] ];
        const Point3D & lB = mVertices[ lTri.mIdxV[1] ];
        const Point3D & lC = mVertices[ lTri.mIdxV[2] ];

        // Compute intersection
        Point3D lI;
		if( pLine.Intersection_NoExc(Plane3D(lA, lB, lC), lI) )
		{
			// Compute local coordinates
			const Vector3D lX = lB - lA;
			const Vector3D lY = lC - lA;
			const Vector3D lU = lI - lA;

			// Fill right hand side
			B.Elt( 0 ) = lX * lU;
			B.Elt( 1 ) = lY * lU;

			// Check if inverse is available
			const Matrix & M_inv = pTriMats[ lIdx ];
			if( M_inv.Rows() == 2 )
			{
				M_inv.Mult( B, XY );

                // Check local coordinates
                // ProstFlex: this code is missing some line-mesh intersections...
//                if( XY.Elt( 0 )>=0 && XY.Elt( 1 )>=0 && (XY.Elt( 0 ) + XY.Elt( 1 ))<=1 )
//					pInters.AddTail( lI );
                // ProstFlex bug: Need to add a tolerance margin to capture all entries
                if( XY.Elt( 0 )>=-1e-6 && XY.Elt( 1 )>=-1e-6 && (XY.Elt( 0 ) + XY.Elt( 1 ))<=1+1e-6 )
                    lSortedInters.AddIntoIncrList( LzServices::Sortable<Point3D,double>(lI, (lI - lLineRoot)*lLineDir), true );

			}
		}
    }

    // ProstFlex bug: ... but now we need to filter out repeating points!
    BrowseList( iSI, lSortedInters )
        pInters.AddTail( lSortedInters.GetAt(iSI).mT );
}

//================================================================================
void Mesh::GetFacesCentroids( vector<Point3D> & pCentroids, vector<double> & pAreas ) const
{
    // Set containers
    pCentroids.resize( mTriangles.size() );
    pAreas.resize( mTriangles.size() );

    // Fill containers
    for( size_t t=0 ; t<mTriangles.size() ; t++ )
    {
        // Triangle and vertices
        const Triangle & lTri = mTriangles[t];
        const Point3D & lA = mVertices[ lTri.mIdxV[0] ];
        const Point3D & lB = mVertices[ lTri.mIdxV[1] ];
        const Point3D & lC = mVertices[ lTri.mIdxV[2] ];

        // Centroid
        pCentroids[t].X() = (lA.X() + lB.X() + lC.X()) / 3;
        pCentroids[t].Y() = (lA.Y() + lB.Y() + lC.Y()) / 3;
        pCentroids[t].Z() = (lA.Z() + lB.Z() + lC.Z()) / 3;

        // Area
        pAreas[t] = ((lB - lA) ^ (lC - lA)).Norm() * 0.5;
    }
}

//================================================================================
void Mesh::GetFacesCentroids( Vector<Point3D> & pCentroids, Vector<double> & pAreas ) const
{
    // Set containers
    pCentroids.Resize( mTriangles.size() );
    pAreas.Resize( mTriangles.size() );

    // Fill containers
    for( size_t t = 0 ; t < mTriangles.size() ; t++ )
    {
        // Triangle and vertices
        const Triangle & lTri = mTriangles[t];
        const Point3D & lA = mVertices[ lTri.mIdxV[0] ];
        const Point3D & lB = mVertices[ lTri.mIdxV[1] ];
        const Point3D & lC = mVertices[ lTri.mIdxV[2] ];

        // Centroid
        pCentroids[t].X() = ( lA.X() + lB.X() + lC.X() ) / 3;
        pCentroids[t].Y() = ( lA.Y() + lB.Y() + lC.Y() ) / 3;
        pCentroids[t].Z() = ( lA.Z() + lB.Z() + lC.Z() ) / 3;

        // Area
        pAreas[t] = ( ( lB - lA ) ^ ( lC - lA ) ).Norm() * 0.5;
    }
}

//================================================================================
Point3D Mesh::WeightedCentroid() const
{
    // Compute face centroids and areas
    Vector<Point3D> lPoints;
    Vector<double> lAreas;
    GetFacesCentroids( lPoints, lAreas );

    // Compute weighted centroid
    return LzMath::ToolBox::WeightedCentroid( lPoints, lAreas );
}

//================================================================================
void Mesh::GetVolumetricCentroids( vector<Point3D> & pCentroids, vector<double> & pVolumes ) const
{
    // Set containers
    pCentroids.resize( mTriangles.size() );
    pVolumes.resize( mTriangles.size() );

    // Pick a reasonable center for tetrahedrons' 4th vertex
    const Point3D lW = LzMath::ToolBox::Centroid( mVertices );

    // Compute each individual point + volume
    for( size_t t=0 ; t<mTriangles.size() ; t++ )
    {
        // Triangle and vertices
        const Triangle & lTri = mTriangles[t];
        const Point3D & lA = mVertices[ lTri.mIdxV[0] ];
        const Point3D & lB = mVertices[ lTri.mIdxV[1] ];
        const Point3D & lC = mVertices[ lTri.mIdxV[2] ];

        // Centroid
        pCentroids[t].X() = ( lA.X() + lB.X() + lC.X() + lW.X() ) / 4.0;
        pCentroids[t].Y() = ( lA.Y() + lB.Y() + lC.Y() + lW.Y() ) / 4.0;
        pCentroids[t].Z() = ( lA.Z() + lB.Z() + lC.Z() + lW.Z() ) / 4.0;

        // Build plane and comp signed tet volume
        const Vector3D lTriNormal = (lC - lA)^(lB - lA); // AC x AB
        const Vector3D lAtoW = lW - lA; // AW

        // Store signed volume
        pVolumes[t] = (lAtoW * lTriNormal) / 6.0;
    }
}

//================================================================================
void Mesh::GetVolumetricCentroids( Vector<Point3D> & pCentroids, Vector<double> & pVolumes ) const
{
    // Set containers
    pCentroids.Resize( mTriangles.size() );
    pVolumes.Resize( mTriangles.size() );

    // Pick a reasonable center for tetrahedrons' 4th vertex
    const Point3D lW = LzMath::ToolBox::Centroid( mVertices );

    // Compute each individual point + volume
    for( size_t t=0 ; t<mTriangles.size() ; t++ )
    {
        // Triangle and vertices
        const Triangle & lTri = mTriangles[t];
        const Point3D & lA = mVertices[ lTri.mIdxV[0] ];
        const Point3D & lB = mVertices[ lTri.mIdxV[1] ];
        const Point3D & lC = mVertices[ lTri.mIdxV[2] ];

        // Centroid
        pCentroids[t].X() = ( lA.X() + lB.X() + lC.X() + lW.X() ) / 4;
        pCentroids[t].Y() = ( lA.Y() + lB.Y() + lC.Y() + lW.Y() ) / 4;
        pCentroids[t].Z() = ( lA.Z() + lB.Z() + lC.Z() + lW.Z() ) / 4;

        // Build plane and comp signed tet volume
        const Vector3D lTriNormal = (lB - lA)^(lC - lA);
        const Vector3D lAtoW = lW - lA;

        // Store signed volume
        pVolumes[t] = (lAtoW * lTriNormal) / 6.0;
    }
}

//================================================================================
Point3D Mesh::VolumetricWeightedCentroid() const
{
    // Compute volumetric centroids and volumes
    Vector<Point3D> lPoints;
    Vector<double> lVolumes;
    GetVolumetricCentroids( lPoints, lVolumes );

    // Compute weighted centroid
    return LzMath::ToolBox::WeightedCentroid( lPoints, lVolumes );
}

//================================================================================
size_t Mesh::GetBestMirrorPlane( const Vector<Plane3D> & pPlanes ) const
{
    // Check
    if( pPlanes.Size() == 0 )
        LzLogException("", "Cannot find best plane among empty planes set!")

	// Build FDC
	FastDistComp lFDC;
    lFDC.Set( *this, LzTriModel::Patching::Split );

	// Find minimal mean error
	double lMinMeanError = std::numeric_limits<double>::max();
    size_t lBestPlane;
    for( size_t p=0 ; p<pPlanes.Size() ; p++ )
	{
		// Get on the plane
		const Plane3D & lP = pPlanes[p];

		// Compute mean error
		double lMeanError = 0.0;
		Point3D lTmp;
		Vector3D lTmp2;
        for( size_t v=0 ; v<mVertices.size() ; v++ )
			lMeanError += fabs( lFDC.SignedDistance(Tree3D::NearestMode::Accurate,
													lP.Symmetrical(mVertices[v]), lTmp, lTmp2) );
//*** DEBUG
//*** DEBUG
{
    LzLogM("", "Mean error for plane "<<lP.ToString_ABCD()<<" = "<<lMeanError);
	Mesh lTmp = *this;
	lTmp.Flip( lP );
    lTmp.Save(LzServices::StartUpPath() + "/_tmp_MirroredMesh_" + std::to_string(p) + ".obj");
}
//*** DEBUG
//*** DEBUG

		// Update
		if( lMinMeanError > lMeanError )
		{
			lMinMeanError = lMeanError;
			lBestPlane = p;
		}
	}

	return lBestPlane;
}

//================================================================================
size_t Mesh::NearestVertexIndex( const Point3D & pPt ) const
{
    // Check
    if( mVertices.size() == 0 )
        LzLogException("", "Cannot find nearest vertex index in an empty mesh!")

    double lMinDist = +std::numeric_limits<double>::max();
    size_t lMinIdx;

    for( size_t v=0 ; v<mVertices.size() ; v++ )
    {
        const double lDist = pPt.DistanceTo( mVertices[v] );

        if( lMinDist > lDist )
        {
            lMinDist = lDist;
            lMinIdx = v;
        }
    }

    return lMinIdx;
}

//================================================================================
size_t Mesh::GetExtremumIndex( const Vector3D & pAlongDir ) const
{
    // Check
    if( mVertices.size() == 0 )
        LzLogException("", "Cannot find extremum vertex index in an empty mesh!")

    double lMaxScal = -std::numeric_limits<double>::max();
    size_t lExtIdx;

    for( size_t v=0 ; v<mVertices.size() ; v++ )
    {
        const double Scal = pAlongDir*(mVertices[v] - Point3D{0,0,0});

        if( lMaxScal < Scal )
        {
            lMaxScal = Scal;
            lExtIdx = v;
        }
    }

    return lExtIdx;

}
#pragma endregion


#pragma region "Properties"
//================================================================================
double Mesh::Volume() const
{
    double lVol = 0;

    Point3D lD = LzMath::ToolBox::Centroid( mVertices );

    for( size_t t = 0 ; t < mTriangles.size() ; t++ )
    {
        const Triangle & lT = mTriangles[t];

        const Point3D & lA = mVertices[ lT.mIdxV[2] ];
        const Point3D & lB = mVertices[ lT.mIdxV[1] ];
        const Point3D & lC = mVertices[ lT.mIdxV[0] ];

        const Vector3D lAB = lB - lA;
        const Vector3D lAC = lC - lA;
        const Vector3D lAD = lD - lA;

        double lTetVol = ( lAB ^ lAC ) * lAD / 6;

        // Update
        lVol += lTetVol;
    }

    return lVol;
}

//================================================================================
double Mesh::Area() const
{
#if 0
    double lSurface = 0.0;

    for( size_t iT = 0 ; iT < pMesh.mTriangles.size(); iT++ )
    {
        const Triangle & lT = pMesh.mTriangles[iT];
        const Point3D  & lPt0 = pMesh.mVertices[ lT.mIdxV[0] ];
        const Point3D  & lPt1 = pMesh.mVertices[ lT.mIdxV[1] ];
        const Point3D  & lPt2 = pMesh.mVertices[ lT.mIdxV[2] ];

        // Compute triangle surface, using Heron formula
        const double lA = lPt0.DistanceTo( lPt1 );
        const double lB = lPt1.DistanceTo( lPt2 );
        const double lC = lPt2.DistanceTo( lPt0 );
        const double lS = 0.5 * ( lA + lB + lC );

        lSurface += sqrt( lS * ( lS - lA ) * ( lS - lB ) * ( lS - lC ) );
    }

    return lSurface;
#else
    double lArea = 0;

    for( size_t t = 0 ; t < mTriangles.size() ; t++ )
    {
        const Triangle & lT = mTriangles[t];

        const Point3D & lA = mVertices[ lT.mIdxV[0] ];
        const Point3D & lB = mVertices[ lT.mIdxV[1] ];
        const Point3D & lC = mVertices[ lT.mIdxV[2] ];

        const Vector3D lAB = lB - lA;
        const Vector3D lAC = lC - lA;

        double lTriArea = (lAB ^ lAC).Norm();

        // Update
        lArea += lTriArea;
    }

    return lArea / 2;
#endif
}

//================================================================================
//bool Mesh::ProbeInversion( DistanceComputer::Patching pPatch ) const
//{
//  // Analyze topology
//  MeshTopology lTopo;
//  lTopo.SetMesh( *this );
//
//  // Check
//  if( lTopo.Cracks().size() )
//      LzLogErr("", "Inversion test might be meaningless: mesh has cracks!");
//
//  // Check
//  if( lTopo.MultiEdges().size() )
//      LzLogErr("", "Inversion test might be meaningless: mesh has multi-edges!");
//
//  // Set distance computer
//  DistanceComputer lDC;
//  lDC.SetTopology( lTopo, pPatch );
//
//  // Compute external probe points
//  BBox lBBox( mVertices );
//  Point3D lExtPts[6] =
//  {
//      lBBox.Center()+Vector3D( lBBox.Size(0), 0, 0 ), lBBox.Center()-Vector3D( lBBox.Size(0), 0, 0 ),
//      lBBox.Center()+Vector3D( lBBox.Size(1), 0, 0 ), lBBox.Center()-Vector3D( lBBox.Size(1), 0, 0 ),
//      lBBox.Center()+Vector3D( lBBox.Size(2), 0, 0 ), lBBox.Center()-Vector3D( lBBox.Size(2), 0, 0 ),
//  };
//
//  // Check distance on external probe points
//  for( int d=0 ; d<6 ; d++ )
//  {
//      Point3D lUnused1;
//      Vector3D lUnused2;
//      if( lDC.SignedDistanceToMesh(lExtPts[d],lUnused1,lUnused2) < 0 )
//          return true;
//  }
//
//  // All clear
//  return false;
//}
#pragma endregion


//    {
//        const Point3D & lV = mVertices[v];
//        glVertex3dv( lV.mV );
//    }
//}
//glEnd();
//
// /#pragma region "Draw"
// //================================================================================
// void Mesh::ResetDisplayList()
// {
// #if !defined(NO_OPENGL)
//     BrowseList( iL, mDLs )
// 		glDeleteLists( mDLs.GetAt(iL), 1 );
// 	mDLs.DelAll();
// #endif
// }

// #if !defined(NO_OPENGL)
// //================================================================================
// void Mesh::Draw()
// {
// #ifndef QT_OPENGL_ES_2 // OpenGL 2.0 draw
//     // Have a display list ?
//     if( mDLs.Count() )
// 	{
// 		BrowseList( iL, mDLs )
// 	        glCallList( mDLs.GetAt(iL) );
// 	}
// 	else
//     {
//         // Has any triangles ?
//         if( mTriangles.empty() )
//         {
// //******* Cette constante devrait etre initialisee en interrogeant l'implementation OpenGL
// static const size_t sMAX_POINTS_IN_LIST = 500000;
// //******* Cette constante devrait etre initialisee en interrogeant l'implementation OpenGL

// 			// Create list by pieces
//             size_t iCountSoFar = 0;
// 			while( iCountSoFar < mVertices.size() )
// 			{
// 				// Alloc new display lists
//                 size_t lDL;
// 				lDL = glGenLists( 1 );
// 				if( !lDL )
// 				{
//                     LzLogErr("",  "Could not alloc a new display list!" );
// 					return;
// 				}

// 				// Feed list
// 				glNewList( lDL, GL_COMPILE_AND_EXECUTE );

// 				// Draw model
// 				glBegin( GL_POINTS );
//                 size_t v = 0;
// 				{
// 					while( v<sMAX_POINTS_IN_LIST && (iCountSoFar + v)<mVertices.size() )
// 					{
// 						const Point3D & lV = mVertices[iCountSoFar + v];
// 						glVertex3dv( lV.mV );

// 						// Next
// 						v++;
// 					}
// 				}
// 				glEnd();

// 				// End list
// 				glEndList();

// 				// Move cursor
// 				iCountSoFar += v;

// 				// Stash DL
// 				mDLs.AddTail( lDL );
// 			}

// //// Alloc new display lists
// //mDL = glGenLists( 1 );
// //if( !mDL )
// //{
// //    LzLogErr("",  "Could not alloc a new display list!" );
// //    return;
// //}
// //
// //// Log
// //LzLogMsg("",  "New DL: " << mDL << "." );
// //
// //// Start list
// //glNewList( mDL, GL_COMPILE_AND_EXECUTE );
// //
// //// Draw model
// //glBegin( GL_POINTS );
// //{
// //    for( size_t v = 0 ; v < mVertices.size() ; v++ )
// /// End list
// //glEndList();
//         }
//         else
//         {
// //******* Cette constante devrait etre initialisee en interrogeant l'implementation OpenGL
// //static const size_t sMAX_TRIANGLES_IN_LIST = 500000;
// static const size_t sMAX_TRIANGLES_IN_LIST = 200000;
// //static const size_t sMAX_TRIANGLES_IN_LIST = 10000;
// //static const size_t sMAX_TRIANGLES_IN_LIST = 2000;
// //******* Cette constante devrait etre initialisee en interrogeant l'implementation OpenGL

// 			// Create list by pieces
//             size_t iCountSoFar = 0;
// 			while( iCountSoFar < mTriangles.size() )
// 			{
// 				// Alloc new display lists
//                 size_t lDL;
// 				lDL = glGenLists( 1 );
// 				if( !lDL )
// 				{
//                     LzLogErr("",  "Could not alloc a new display list!" );
// 					return;
// 				}

// 				// Log
//                 LzLogMsg("",  "New DL: " << lDL << "." );

// 				// Feed list
//                 glNewList( lDL, GL_COMPILE_AND_EXECUTE );

// 				// Draw model
//                 size_t t = 0;
// 				glBegin( GL_TRIANGLES );
// 				{
// 					while( t<sMAX_TRIANGLES_IN_LIST && (iCountSoFar + t)<mTriangles.size() )
// 					{
// 						const Triangle & lTri = mTriangles[iCountSoFar + t];

// 						for( int i=0 ; i<3 ; i++ )
// 						{
// 							// Have a tex coord?
// 							if( mTexCoords2D.Size() )
// 							{
// 								const TexCoord2D & lST = mTexCoords2D[ lTri.mIdxT[i] ];
// 								glTexCoord2d( lST.mS, lST.mT );
// 							}

// 							// Normal and vertex
// 							glNormal3dv( mNormals[lTri.mIdxN[i]].mV );
// 							glVertex3dv( mVertices[lTri.mIdxV[i]].mV );
// 						}

// 						// Next
// 						t++;
// 					}
// 				}
// 				glEnd();

// 				// End list
// 				glEndList();

// 				// Move cursor
// 				iCountSoFar += t;

// 				// Stash DL
// 				mDLs.AddTail( lDL );
// 			}
// /*
// //// Get max number of triangles in DL
// //static const size_t sMaxTrianglesInDL = 900000;
// ////===> no such thing as glGetIntegerv( GL_MAX_TRIANGLES_IN_LIST, &lMaxTrianglesInDL );
// //
// //// Actual number of triangles
// //const size_t lTriCount = mTriangles.size();
// //
// //// Can we turn it into a display list?
// //if( lTriCount < sMaxTrianglesInDL )
// //{
// //    // Alloc new display lists
// //    mDL = glGenLists( 1 );
// //    if( !mDL )
// //    {
// //        LzLogErr("",  "Could not alloc a new display list!" );
// //        return;
// //    }
// //
// //    // Log
// //    LzLogMsg("",  "New DL: " << mDL << "." );
// //
// //    // Start list
// //    glNewList( mDL, GL_COMPILE_AND_EXECUTE );
// //}
// //
// //// Draw model
// //glBegin( GL_TRIANGLES );
// //{
// //    for( size_t t = 0 ; t < lTriCount ; t++ )
// //    {
// //        const Triangle & lTri = mTriangles[t];
// //
// //        for( int i = 0 ; i < 3 ; i++ )
// //        {
// //            // Have a tex coord?
// //            if( mTexCoords2D.Size() )
// //            {
// //                const TexCoord2D & lST = mTexCoords2D[ lTri.mIdxT[i] ];
// //                glTexCoord2d( lST.mS, lST.mT );
// //            }
// //
// //            // Normal and vertex
// //            glNormal3dv( mNormals[lTri.mIdxN[i]].mV );
// //            glVertex3dv( mVertices[lTri.mIdxV[i]].mV );
// //        }
// //    }
// //}
// //glEnd();
// //
// //if( lTriCount < sMaxTrianglesInDL )
// //{
// //    // End list
// //    glEndList();
// //}*/
//         }

//         // Final check
//         GLenum lErrEnum = glGetError();
//         if( lErrEnum != GL_NO_ERROR )
//         {
//             // Release display list
//             ResetDisplayList();

//             // Log
//             LzLogErr("", "OpenGL error number: "<<lErrEnum)
//             LzLogErr("", "OpenGL error string: "<<(const char*)gluErrorString(lErrEnum))
//             LzLogException("", "Too many triangles in display list!" ) // actual error: graphics board out of memory
//         }
//     }
// #else // OpenGL ES 2 Draw

//     LzLogException("", "*** TODO  Mesh::Draw()");

// #endif
// }
// #endif
// #pragma endregion
}
