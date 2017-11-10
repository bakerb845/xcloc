#define _BSD_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <libgen.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/tree.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>
#include <libxml/parser.h>
#include "xcloc_xdmf.h"

#define XML_ENCODING "utf-8"
#define XDMF_VERSION "3"
#define TOPOLOGY3D_TYPE "3DCoRectMesh"
#define TOPOLOGY2D_TYPE "2DCoRectMesh"
#define GEOMETRY3D_TYPE "ORIGIN_DXDYDZ"
#define GEOMETRY2D_TYPE "ORIGIN_DXDY"
#define CHUNK 256

static int writeH5DataItem3d(const int nz, const int ny, const int nx, 
                             const enum xclocPrecision_enum precision,
                             const char *h5File, const char *dataSet,
                             xmlTextWriterPtr *writer);
static int writeGridSpacing3d(const double dz, const double dy, const double dx, 
                              xmlTextWriterPtr *writer);
static int writeOrigin3d(const double z0, const double y0, const double x0,
                         xmlTextWriterPtr *writer);
static int writeTopology3d(const int nz, const int ny, const int nx,
                           xmlTextWriterPtr *writer);

#define CHECKRC(rc) \
 if (rc < 0) \
 { \
     fprintf(stderr, "%s: XML failure on line %d\n", __func__, __LINE__); \
     return -1; \
 }

/*!
 * @brief Frees the XDMF grid strutcure.
 * 
 * @param[out] xdmf   On exit the memory on the xdmf grid structure has
 *                    been released and variables set to 0.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under MIT.
 *
 */
int xcloc_xdmfGrid_free(struct xclocXDMFGrid_struct *xdmf)
{
    int i;
    if (!xdmf->linit){return 0;}
    if (xdmf->h5flNames != NULL)
    {
        for (i=0; i<xdmf->nDataSets; i++)
        {
            if (xdmf->h5flNames[i] != NULL){free(xdmf->h5flNames[i]);}
        }
        free(xdmf->h5flNames);
    }
    if (xdmf->dataSets != NULL)
    {
        for (i=0; i<xdmf->nDataSets; i++)
        {
            if (xdmf->dataSets[i] != NULL){free(xdmf->dataSets[i]);}
        }
        free(xdmf->dataSets);
    }
    if (xdmf->precision){free(xdmf->precision);}
    memset(xdmf, 0, sizeof(struct xclocXDMFGrid_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Initializes the XDMF structure.  This structure will collect
 *        like gridded models.
 *
 * @param[in] nx      Number of x grid points in model.  In 3D this is the
 *                    fastest changing model dimension.
 * @param[in] ny      Number of y grid points in model. 
 * @param[in] nz      Number of z grid points in model.  In 3D this is the
 *                    slowest chaning model dimension.
 *                    varying model dimension.
 * @param[in] dx      Grid spacing (meters) in x.
 * @param[in] dy      Grid spacing (meters) in y.
 * @param[in] dz      Grid spacing (meters) in z.
 * @param[in] x0      x model origin (meters).
 * @param[in] y0      y model origin (meters). 
 * @param[in] z0      z model origin (meters).
 *
 * @param[out] xdmf   Initialized XDMF gridded model structure.
 *
 * @copyright Ben Baker distributed under MIT.
 *
 */
int xcloc_xdmfGrid_initialize(
    const int nx, const int ny, const int nz,
    const double dx, const double dy, const double dz,
    const double x0, const double y0, const double z0,
    struct xclocXDMFGrid_struct *xdmf)
{
    size_t sizeEnum = sizeof(enum xclocPrecision_enum);
    memset(xdmf, 0, sizeof(struct xclocXDMFGrid_struct));
    if (nx < 1 || ny < 1 || nz < 1 ||
        (nx > 1 && dx == 0.0) || (ny > 1 && dy == 0.0) || (nz > 1 && dz == 0.0))
    {
        if (nx < 1){fprintf(stderr, "%s: nx must be positive\n", __func__);}
        if (ny < 1){fprintf(stderr, "%s: ny must be positive\n", __func__);}
        if (nz < 1){fprintf(stderr, "%s: nz must be positive\n", __func__);}
        if (nx > 1 && dx == 0.0)
        {
            fprintf(stderr, "%s: dx cannot be 0\n", __func__);
        }
        if (ny > 1 && dy == 0.0)
        {
            fprintf(stderr, "%s: dy cannot be 0\n", __func__);
        }
        if (nz > 1 && dz == 0.0)
        {
            fprintf(stderr, "%s: dz cannot be 0\n", __func__);
        }
        return -1;
    }
    if ((nx == 1 && ny == 1) || (nx == 1 && nz == 1) || (ny == 1 && nz == 1))
    {
        fprintf(stderr, "%s: Can't deal with 1d model\n", __func__);
        return -1;
    }
    if (nx == 1 || ny == 1 || nz == 1)
    {
        xdmf->lis2d = true;
        if (nx == 1)
        {
            xdmf->nx = ny;
            xdmf->ny = nz;
            xdmf->nz = 1;
            xdmf->x0 = y0;
            xdmf->y0 = z0;
            xdmf->z0 = 0.0;
            xdmf->dx = dy;
            xdmf->dy = dz;
            xdmf->dz = 1.0;
        }
        else if (ny == 1)
        {
            xdmf->nx = nx;
            xdmf->ny = nz;
            xdmf->nz = 1;
            xdmf->x0 = x0;
            xdmf->y0 = z0;
            xdmf->z0 = 0.0;
            xdmf->dx = dx;
            xdmf->dy = dz;
            xdmf->dz = 1.0;
        }
        else
        {
            xdmf->nx = nx;
            xdmf->ny = ny;
            xdmf->nz = 1;
            xdmf->x0 = x0;
            xdmf->y0 = y0;
            xdmf->z0 = 0.0;
            xdmf->dx = dx;
            xdmf->dy = dy;
            xdmf->dz = 1.0;
        }
    }
    else
    {
        xdmf->nx = nx;
        xdmf->ny = ny;
        xdmf->nz = nz;
        xdmf->x0 = x0;
        xdmf->y0 = y0;
        xdmf->z0 = z0;
        xdmf->dx = dx;
        xdmf->dy = dy;
        xdmf->dz = dz;
    }
    xdmf->mDataSets = CHUNK;
    xdmf->h5flNames = (char **) calloc(xdmf->mDataSets, sizeof(char *));
    xdmf->dataSets  = (char **) calloc(xdmf->mDataSets, sizeof(char *));
    xdmf->precision = (enum xclocPrecision_enum *)
                      calloc(xdmf->mDataSets, sizeEnum);
    xdmf->linit = true;
    return 0;
}
//============================================================================//
/*!
 * @brief Appends the HDF5 archive file and dataset  in the archive file to 
 *        the XDMF structure.
 *
 * @param[in] h5flName   Name of HDF5 file.
 * @param[in] dataSet    Absolute path to the dataset in the HDF5 file archive.
 * @param[in] precision  Precision of dataset.
 *
 * @param[in,out] xdmf   On input this is the initialized grid. \n
 *                       On exit the dataset in the HDF5 file has been 
 *                       appended.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xdmfGrid_add(const char *h5flName, const char *dataSet,
                       const enum xclocPrecision_enum precision,
                       struct xclocXDMFGrid_struct *xdmf)
{
    char **work1, **work2;
    char *h5b, *temp;
    enum xclocPrecision_enum *prec;
    int i;
    temp = strdup(h5flName);
    h5b = basename(temp);
    if (xdmf->nDataSets + 1 > xdmf->mDataSets)
    {
        work1 = (char **) calloc(xdmf->mDataSets + CHUNK, sizeof(char *));
        work2 = (char **) calloc(xdmf->mDataSets + CHUNK, sizeof(char *));
        prec = (enum xclocPrecision_enum *)
               calloc(xdmf->mDataSets + CHUNK,
                      sizeof(enum xclocPrecision_enum));
        for (i=0; i<xdmf->nDataSets; i++)
        {
            prec[i] = xdmf->precision[i];
            work1[i] = strdup(xdmf->h5flNames[i]);
            work2[i] = strdup(xdmf->dataSets[i]);
            free(xdmf->h5flNames[i]);
            free(xdmf->dataSets[i]);
        }
        free(xdmf->precision);
        free(xdmf->h5flNames);
        free(xdmf->dataSets);
        xdmf->precision = prec;
        xdmf->h5flNames = work1;
        xdmf->dataSets  = work2;
        xdmf->mDataSets = xdmf->mDataSets + CHUNK;
    }
    i = xdmf->nDataSets;
    xdmf->precision[i] = precision;
    xdmf->h5flNames[i] = strdup(h5b);
    xdmf->dataSets[i]  = strdup(dataSet);
    xdmf->nDataSets = xdmf->nDataSets + 1;
    free(temp);
    return 0; 
}
//============================================================================//
/*!
 * @brief Writes the XDMF grid file.
 *
 * @param[in] xdmfFile   Name of XDMF file to write.  This must exist
 *                       alongside the HDF5 archive file(s).
 * @param[in] xdmf       The initialized XDMF grid structure.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xdmf_writeGrid(const char *xdmfFile,
                         const struct xclocXDMFGrid_struct xdmf)
{
    FILE *xmlFileHandle;
    xmlTextWriterPtr writer;
    xmlBufferPtr buf;
    char *xmlmsg;
    size_t msglen;
    int i, rc;
    if (!xdmf.linit)
    {
        fprintf(stderr, "%s: xdmf structure not initialized\n", __func__);
        return -1;
    }
    if (xdmf.nDataSets == 0)
    {
        fprintf(stdout, "%s: No datasets - leaving early\n", __func__);
        return 0;
    }
    // Create a new XML buffer to which the XML document will be written
    buf = xmlBufferCreate();
    if (buf == NULL)
    {
        fprintf(stderr, "%s: Error creating XML buffer!", __func__);
        return -1;
    }
    // Create a new xmlWriter for uri with no compression
    writer = xmlNewTextWriterMemory(buf, 0);
    if (writer == NULL)
    {
        fprintf(stderr, "%s: Error creating xml writer", __func__);
        return -1;
    }
    xmlTextWriterSetIndentString(writer, "  "); //"\t");
    // Start the document with default xml version
    rc = xmlTextWriterStartDocument(writer, NULL, XML_ENCODING, NULL);
    if (rc < 0)
    {
        fprintf(stderr, "%s: Error starting writer", __func__);
        return -1;
    }
    // <Xdmf>
    xmlTextWriterSetIndent(writer, 1);
    rc = xmlTextWriterStartElement(writer, BAD_CAST "Xdmf\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "xmlns:xi\0",
                                     BAD_CAST "http://www.w3.org/2003/XInclude\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Version\0",
                                     BAD_CAST XDMF_VERSION);
    // <Domain>
    rc = xmlTextWriterStartElement(writer, BAD_CAST "Domain\0"); CHECKRC(rc);
    // <Grid>
    if (xdmf.nDataSets > 0)
    {
        // <Grid>
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Grid\0"); CHECKRC(rc);
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name\0",
                                         BAD_CAST "DSM Images\0"); CHECKRC(rc);
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "GridType\0",
                                         BAD_CAST "Uniform\0"); CHECKRC(rc);
        // <Topology>
        if (xdmf.lis2d) 
        {
            writeTopology3d(1, xdmf.ny, xdmf.nx, &writer);
        }
        else 
        {
            writeTopology3d(xdmf.nz, xdmf.ny, xdmf.nx, &writer);
        }
        // </Topology>
        // <Geometry>
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Geometry\0");
        CHECKRC(rc);
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "GeometryType\0",
                                         BAD_CAST GEOMETRY3D_TYPE); CHECKRC(rc);
        if (xdmf.lis2d)
        {
            writeOrigin3d(0.0, xdmf.y0, xdmf.x0, &writer);
            writeGridSpacing3d(0.0, xdmf.dy, xdmf.dx, &writer);
        }
        else
        {
            writeOrigin3d(xdmf.z0, xdmf.y0, xdmf.x0, &writer);
            writeGridSpacing3d(xdmf.dz, xdmf.dy, xdmf.dx, &writer);
        }
        //  </Geometry>
        rc = xmlTextWriterEndElement(writer); CHECKRC(rc);
    }
    // Write the data sets in the collection
    for (i=0; i<xdmf.nDataSets; i++)
    {
        //  <Attribute>
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Attribute\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name\0",
                                         BAD_CAST xdmf.dataSets[i]);
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Type\0",
                                         BAD_CAST "Scalar\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Center\0",
                                         BAD_CAST "Node\0");
        //    <DataItem>
        if (xdmf.lis2d)
        {
            rc = writeH5DataItem3d(1, xdmf.ny, xdmf.nx, xdmf.precision[i],
                                   xdmf.h5flNames[i], xdmf.dataSets[i],
                                   &writer); CHECKRC(rc);
        }
        else
        {
            rc = writeH5DataItem3d(xdmf.nz, xdmf.ny, xdmf.nx, xdmf.precision[i],
                                   xdmf.h5flNames[i], xdmf.dataSets[i],
                                   &writer); CHECKRC(rc);
        } //  </DataItem>
        //  </Attribute>
        rc = xmlTextWriterEndElement(writer); CHECKRC(rc);
    } // Loop on datasets
    // </Grid>
    if (xdmf.nDataSets > 1){rc = xmlTextWriterEndElement(writer); CHECKRC(rc);}
    rc = xmlTextWriterEndElement(writer); CHECKRC(rc); //</Domain>
    rc = xmlTextWriterEndElement(writer); CHECKRC(rc); //</Xdmf>
    // Finally copy the char * XML message
    xmlFreeTextWriter(writer);
    xmlCleanupCharEncodingHandlers();
    msglen = (size_t) xmlStrlen(buf->content);
    xmlmsg = (char *) calloc(msglen + 1, sizeof(char));
    strncpy(xmlmsg, (const char *)buf->content, msglen);
    // Write it 
    xmlFileHandle = fopen(xdmfFile, "w");
    fprintf(xmlFileHandle, xmlmsg);
    fclose(xmlFileHandle);
    // Clean up
    xmlCleanupParser();
    xmlBufferFree(buf);
    xmlDictCleanup();
    xmlCleanupThreads();
    free(xmlmsg);
    return 0;
}
//============================================================================//

/*!
 * @brief Convenience function to write the origin in 3d.
 */
static int writeOrigin3d(const double z0, const double y0, const double x0,
                         xmlTextWriterPtr *writer)
{
    char cori[128];
    int rc;
    // Origin
    rc = xmlTextWriterStartElement(*writer, BAD_CAST "DataItem\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Name\0",
                                     BAD_CAST "Origin\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST "3\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "DataType\0",
                                     BAD_CAST "Float\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Precision\0",
                                     BAD_CAST "8\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Format\0",
                                     BAD_CAST "XML\0");
    memset(cori, 0, 128*sizeof(char));
    sprintf(cori, "%e %e %e", z0, y0, x0);
    rc = xmlTextWriterWriteString(*writer, BAD_CAST cori);
    rc = xmlTextWriterEndElement(*writer); //</DataItem>
    return rc;
}
/*
static int writeOrigin2d(const double y0, const double x0,
                         xmlTextWriterPtr *writer)
{
    char cori[128];
    int rc;
    // Origin
    rc = xmlTextWriterStartElement(*writer, BAD_CAST "DataItem\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Name\0",
                                     BAD_CAST "Origin\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST "2\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "DataType\0",
                                     BAD_CAST "Float\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Precision\0",
                                     BAD_CAST "8\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Format\0",
                                     BAD_CAST "XML\0");
    memset(cori, 0, 128*sizeof(char));
    sprintf(cori, "%e %e", y0, x0);
    rc = xmlTextWriterWriteString(*writer, BAD_CAST cori);
    rc = xmlTextWriterEndElement(*writer); //</DataItem>
    return rc;
}
*/
//============================================================================//
/*!
 * @brief Convenience function to write the topology in 3d.
 */
static int writeTopology3d(const int nz, const int ny, const int nx,
                           xmlTextWriterPtr *writer)
{
    char cdim[128];
    int rc;
    memset(cdim, 0, 128*sizeof(char));
    sprintf(cdim, "%d %d %d", nz, ny, nx);
    // Topology
    rc = xmlTextWriterStartElement(*writer, BAD_CAST "Topology\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "TopologyType\0",
                                     BAD_CAST TOPOLOGY3D_TYPE);
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST cdim);
    rc = xmlTextWriterEndElement(*writer); //</Topology>
    return rc;
}
//============================================================================//
/*
static int writeTopology2d(const int ny, const int nx,
                         xmlTextWriterPtr *writer)
{
    char cdim[128];
    int rc;
    memset(cdim, 0, 128*sizeof(char));
    sprintf(cdim, "%d %d", ny, nx);
    // Topology
    rc = xmlTextWriterStartElement(*writer, BAD_CAST "Topology\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "TopologyType\0",
                                     BAD_CAST TOPOLOGY2D_TYPE);
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST cdim);
    rc = xmlTextWriterEndElement(*writer); //</Topology>
    return rc;
}
*/
//============================================================================//
/*!
 * @brief Convenience function to write the grid spacing in 3d.
 */
static int writeGridSpacing3d(const double dz, const double dy, const double dx,
                              xmlTextWriterPtr *writer)
{
    char cspace[128];
    int rc;
    memset(cspace, 0, 128*sizeof(char));
    sprintf(cspace, "%e %e %e", dz, dy, dx);
    // Grid spacing
    rc = xmlTextWriterStartElement(*writer, BAD_CAST "DataItem\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Name\0",
                                     BAD_CAST "Spacing\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST "3\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "DataType\0",
                                     BAD_CAST "Float\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Precision\0",
                                     BAD_CAST "8\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Format\0",
                                     BAD_CAST "XML\0");
    rc = xmlTextWriterWriteString(*writer, BAD_CAST cspace);
    rc = xmlTextWriterEndElement(*writer); //</DataItem>
    return rc;
}
//============================================================================//
/*
static int writeGridSpacing2d(const double dy, const double dx, 
                              xmlTextWriterPtr *writer)
{
    char cspace[128];
    int rc; 
    memset(cspace, 0, 128*sizeof(char));
    sprintf(cspace, "%e %e", dy, dx);
    // Grid spacing
    rc = xmlTextWriterStartElement(*writer, BAD_CAST "DataItem\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Name\0",
                                     BAD_CAST "Spacing\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST "2\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "DataType\0",
                                     BAD_CAST "Float\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Precision\0",
                                     BAD_CAST "8\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Format\0",
                                     BAD_CAST "XML\0");
    rc = xmlTextWriterWriteString(*writer, BAD_CAST cspace);
    rc = xmlTextWriterEndElement(*writer); //</DataItem>
    return rc; 
}
*/
//============================================================================//

static int writeH5DataItem3d(const int nz, const int ny, const int nx,
                             const enum xclocPrecision_enum precision,
                             const char *h5File, const char *dataSet,
                             xmlTextWriterPtr *writer)
{
    char path[PATH_MAX];
    char cdim[128];
    int rc;
    // DataItem
    memset(path, 0, PATH_MAX*sizeof(char));
    memset(cdim, 0, 128*sizeof(char));
    sprintf(cdim, "%d %d %d", nz, ny, nx);
    rc = xmlTextWriterStartElement(*writer, BAD_CAST "DataItem\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "DataType\0",
                                     BAD_CAST "Float\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST cdim);
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Format\0",
                                     BAD_CAST "HDF\0");
    if (precision == XCLOC_SINGLE_PRECISION)
    {
        rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Precision\0",
                                         BAD_CAST "4\0");
    }
    else
    {
        rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Precision\0",
                                         BAD_CAST "8\0");
    }
    sprintf(path, "%s:%s", h5File, dataSet);
    rc = xmlTextWriterWriteString(*writer, BAD_CAST path);
    //   </DataItem>
    rc = xmlTextWriterEndElement(*writer); //</DataItem>
    return rc;
}
//============================================================================//
/*
static int writeH5DataItem2d(const int ny, const int nx, 
                             const enum xclocPrecision_enum precision,
                             const char *h5File, const char *dataSet,
                             xmlTextWriterPtr *writer)
{
    char path[PATH_MAX];
    char cdim[128];
    int rc; 
    // DataItem
    memset(path, 0, PATH_MAX*sizeof(char));
    memset(cdim, 0, 128*sizeof(char));
    sprintf(cdim, "%d %d", ny, nx);
    rc = xmlTextWriterStartElement(*writer, BAD_CAST "DataItem\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "DataType\0",
                                     BAD_CAST "Float\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST cdim);
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Format\0",
                                     BAD_CAST "HDF\0");
    if (precision == XCLOC_SINGLE_PRECISION)
    {   
        rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Precision\0",
                                         BAD_CAST "4\0");
    }   
    else
    {   
        rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Precision\0",
                                         BAD_CAST "8\0");
    }   
    sprintf(path, "%s:%s", h5File, dataSet);
    rc = xmlTextWriterWriteString(*writer, BAD_CAST path);
    //   </DataItem>
    rc = xmlTextWriterEndElement(*writer); //</DataItem>
    return rc; 
}
*/
