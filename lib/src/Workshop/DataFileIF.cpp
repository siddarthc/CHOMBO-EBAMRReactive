#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DataFileIF.H"

#include "NamespaceHeader.H"

using std::cin;

DataFileIF::DataFileIF(const DataFileIF::DataType& a_dataType,
                       const Real&                 a_value,
                       const bool&                 a_inside)
{
  // Read an entire header from "cin" - see .H file
  ReadFullHeader(m_num,m_spacing,m_origin,cin);

  // Read all the data from "cin"
  ReadData(m_noDataValue,cin,a_dataType,m_num);

  // Save some state
  m_value = a_value;
  m_inside = a_inside;

  // Make the unit cube IntVectSet
  MakeCorners();
}

DataFileIF::DataFileIF(const char* const           a_filename,
                       const DataFileIF::DataType& a_dataType,
                       const Real&                 a_value,
                       const bool&                 a_inside)
{
  ifstream curFile;

  // Open the named file for input
  OpenFile(curFile,a_filename);

  // Read an entire header from the file - see .H file
  ReadFullHeader(m_num,m_spacing,m_origin,curFile);

  // Read all the data from the file
  ReadData(m_noDataValue,curFile,a_dataType,m_num);

  // Close the named file
  CloseFile(curFile);

  // Save some state
  m_value = a_value;
  m_inside = a_inside;

  // Make the unit cube IntVectSet
  MakeCorners();
}

DataFileIF::DataFileIF(const DataFileIF::DataType& a_dataType,
                       const RealVect&             a_spacing,
                       const RealVect&             a_origin,
                       const Real&                 a_value,
                       const bool&                 a_inside)
{
  // Read an dimension of the data from "cin" - see .H file
  ReadMinHeader(m_num,cin);

  // Save the other "header information"
  m_spacing = a_spacing;
  m_origin = a_origin;

  // Read all the data from "cin"
  ReadData(m_noDataValue,cin,a_dataType,m_num);

  // Save some state
  m_value = a_value;
  m_inside = a_inside;

  // Make the unit cube IntVectSet
  MakeCorners();
}

DataFileIF::DataFileIF(const char* const           a_filename,
                       const DataFileIF::DataType& a_dataType,
                       const RealVect&             a_spacing,
                       const RealVect&             a_origin,
                       const Real&                 a_value,
                       const bool&                 a_inside)
{
  ifstream curFile;

  // Open the named file for input
  OpenFile(curFile,a_filename);

  // Read an dimension of the data from the file - see .H file
  ReadMinHeader(m_num,curFile);

  // Save the other "header information"
  m_spacing = a_spacing;
  m_origin = a_origin;

  // Read all the data from the file
  ReadData(m_noDataValue,curFile,a_dataType,m_num);

  // Close the named file
  CloseFile(curFile);

  // Save some state
  m_value = a_value;
  m_inside = a_inside;

  // Make the unit cube IntVectSet
  MakeCorners();
}

DataFileIF::DataFileIF(const DataFileIF::DataType& a_dataType,
                       const IntVect&              a_num,
                       const RealVect&             a_spacing,
                       const RealVect&             a_origin,
                       const Real&                 a_value,
                       const bool&                 a_inside)
{
  // Save the "header information"
  m_num = a_num;
  m_spacing = a_spacing;
  m_origin = a_origin;

  // Read all the data from "cin"
  ReadData(m_noDataValue,cin,a_dataType,m_num);

  // Save some state
  m_value = a_value;
  m_inside = a_inside;

  // Make the unit cube IntVectSet
  MakeCorners();
}

DataFileIF::DataFileIF(const char* const           a_filename,
                       const DataFileIF::DataType& a_dataType,
                       const IntVect&              a_num,
                       const RealVect&             a_spacing,
                       const RealVect&             a_origin,
                       const Real&                 a_value,
                       const bool&                 a_inside)
{
  ifstream curFile;

  // Open the named file for input
  OpenFile(curFile,a_filename);

  // Save the "header information"
  m_num = a_num;
  m_spacing = a_spacing;
  m_origin = a_origin;

  // Read all the data from the file
  ReadData(m_noDataValue,curFile,a_dataType,m_num);

  // Close the named file
  CloseFile(curFile);

  // Save some state
  m_value = a_value;
  m_inside = a_inside;

  // Make the unit cube IntVectSet
  MakeCorners();
}

DataFileIF::DataFileIF(const DataFileIF& a_inputIF)
{
  // Save a refcounted pointer to the data
  m_ascii_data  = a_inputIF.m_ascii_data;
  m_binary_data = a_inputIF.m_binary_data;
  m_noDataValue = a_inputIF.m_noDataValue;

  // Copy all the other data
  m_num = a_inputIF.m_num;

  m_spacing = a_inputIF.m_spacing;
  m_origin = a_inputIF.m_origin;

  m_value = a_inputIF.m_value;
  m_inside = a_inputIF.m_inside;

  // Make the unit cube IntVectSet
  MakeCorners();
}

DataFileIF::DataFileIF(const RefCountedPtr<FArrayBox>               a_ascii_data,
                       const RefCountedPtr<BaseFab<unsigned char> > a_binary_data,
                       const Real&                                  a_noDataValue,
                       const IntVect&                               a_num,
                       const RealVect&                              a_spacing,
                       const RealVect&                              a_origin,
                       const Real&                                  a_value,
                       const bool&                                  a_inside)
{
  // Save a refcounted pointer to the data
  m_ascii_data  = a_ascii_data;
  m_binary_data = a_binary_data;
  m_noDataValue = a_noDataValue;

  // Copy all the other data
  m_num = a_num;

  m_spacing = a_spacing;
  m_origin = a_origin;

  m_value = a_value;
  m_inside = a_inside;

  // Make the unit cube IntVectSet
  MakeCorners();
}

DataFileIF::~DataFileIF()
{
}

void DataFileIF::GetHeader(IntVect&  a_num,
                           RealVect& a_spacing,
                           RealVect& a_origin) const
{
  // Copy header information over
  a_num = m_num;
  a_spacing = m_spacing;
  a_origin = m_origin;
}

void DataFileIF::GetParams(Real& a_value,
                           bool& a_inside) const
{
  // Copy parameter information over
  a_value = m_value;
  a_inside = m_inside;
}

void DataFileIF::SetParams(const Real& a_value,
                           const bool& a_inside)
{
  // Set parameter information
  m_value = a_value;
  m_inside = a_inside;
}

void DataFileIF::SetNoDataValue(const Real& a_value)
{
  // value to use when we are outside
  m_noDataValue = a_value;
}

Real DataFileIF::value(const RealVect& a_point) const
{
  Real retval;

  // The box of the stored data
  Box dataBox;
  if (m_ascii_data != NULL)
  {
    dataBox = m_ascii_data->box();
  }
  else
  {
    dataBox = m_binary_data->box();
  }

  // The index in the data corresponding to a_point
  IntVect loCorner;
  // The fraction of a cell from the low index to a_point
  RealVect fraction;

  // Determine the low index and the fraction
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    Real index;

    // The floating point "index" of a_point in the data
    index = (a_point[idir] - m_origin[idir]) / m_spacing[idir];

    // The integer portion of that "index"
    loCorner[idir] = (int)floor(index);

    // Make sure the low corner in one inside the high end of the data
    if (loCorner[idir] == dataBox.bigEnd(idir))
    {
      loCorner[idir] -= 1;
    }

    // The fraction portion of the "index"
    fraction[idir] = index - loCorner[idir];
  }

  IntVect hiCorner(loCorner);
  hiCorner += IntVect::Unit;

  // Make a box where data is needed to compute a value
  Box cell(loCorner,hiCorner);

  // If so, use trilinear interpolation to get a return value
  retval = 0.0;

  // Iterate over all the corners
  int comp = 0;
  IVSIterator ivsit(m_corners);
  for (ivsit.begin(); ivsit.ok(); ++ivsit)
  {
    IntVect curIV(loCorner);
    IntVect incIV(ivsit());
    Real curValue;

    // Index of the current data
    curIV += incIV;

    // Check to see if current index in insid the data
    if (dataBox.contains(curIV))
    {
      // If so, get the current data value
      if (m_ascii_data != NULL)
      {
        curValue = (*m_ascii_data)(curIV,comp);
      }
      else
      {
        curValue = (*m_binary_data)(curIV,comp);
      }
    }
    else
    {
      // If not, use the no data value
      // pout() << "Using no data value: " << m_noDataValue << endl;
      curValue = m_noDataValue;
    }

    // Weight it appropriate based on the index fractions computed above
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (incIV[idir] == 0)
      {
        curValue *= 1 - fraction[idir];
      }
      else
      {
        curValue *= fraction[idir];
      }
    }

    // Add it into the total
    retval += curValue;
  }

  // Adjust the level set from zero to m_value
  retval -= m_value;

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* DataFileIF::newImplicitFunction() const
{
  DataFileIF* dataFilePtr = new DataFileIF(m_ascii_data,
                                           m_binary_data,
                                           m_noDataValue,
                                           m_num,
                                           m_spacing,
                                           m_origin,
                                           m_value,
                                           m_inside);

  return static_cast<BaseIF*>(dataFilePtr);
}

void DataFileIF::OpenFile(ifstream&         a_file,
                          const char* const a_filename)
{
  a_file.open(a_filename);

  if (!a_file.good())
  {
    MayDay::Abort("DataFileIF::OpenFile - Unable to open data file");
  }
}

void DataFileIF::CloseFile(ifstream& a_file)
{
  a_file.close();
}

void DataFileIF::ReadMinHeader(IntVect& a_num,
                               istream& a_file)
{
  char curLine[1024];
  int length = 1024;

  // Read first line - number of grid points in each direction
  a_file.getline(curLine,length);

  int num[3];
  int numFound;

  numFound = sscanf(curLine,"%d%d%d",&num[0],&num[1],&num[2]);

  if (numFound < SpaceDim)
  {
    MayDay::Abort("DataFileIF::ReadMinHeader - Unable to read number of grid points in each direction from data file");
  }
  else
  {
    if (numFound > SpaceDim)
    {
      MayDay::Warning("Found more entries in header line than the number of dimensions - wrong dimension data?");
    }

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_num[idir] = num[idir];
    }
  }
}

void DataFileIF::ReadFullHeader(IntVect& a_num,
                                RealVect& a_spacing,
                                RealVect& a_origin,
                                istream&  a_file)
{
  char curLine[1024];
  int length = 1024;

  // Read first line - number of grid points in each direction
  a_file.getline(curLine,length);

  int num[3];
  int numFound;

  numFound = sscanf(curLine,"%d%d%d",&num[0],&num[1],&num[2]);

  if (numFound < SpaceDim)
  {
    MayDay::Abort("DataFileIF::ReadFullHeader - Unable to read number of grid points in each direction from data file");
  }
  else
  {
    if (numFound > SpaceDim)
    {
      MayDay::Warning("Found more entries in header line than the number of dimensions - wrong dimension data?");
    }

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_num[idir] = num[idir];
    }
  }

  // Read second line - grid spacing
  a_file.getline(curLine,length);

  double spacing[3];

  numFound = sscanf(curLine,"%lf%lf%lf",&spacing[0],&spacing[1],&spacing[2]);

  if (numFound < SpaceDim)
  {
    MayDay::Abort("DataFileIF::ReadFullHeader - Unable to read number of grid spacing from data file");
  }
  else
  {
    if (numFound > SpaceDim)
    {
      MayDay::Warning("Found more entries in header line than the number of dimensions - wrong dimension data?");
    }

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_spacing[idir] = spacing[idir];
    }
  }

  // Read third line - grid origin
  a_file.getline(curLine,length);

  double origin[3];

  numFound = sscanf(curLine,"%lf%lf%lf",&origin[0],&origin[1],&origin[2]);

  if (numFound < SpaceDim)
  {
    MayDay::Abort("DataFileIF::ReadFullHeader - Unable to read number of grid origin from data file");
  }
  else
  {
    if (numFound > SpaceDim)
    {
      MayDay::Warning("Found more entries in header line than the number of dimensions - wrong dimension data?");
    }

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_origin[idir] = origin[idir];
    }
  }
}

void DataFileIF::ReadData(Real&                       a_maxValue,
                          istream&                    a_file,
                          const DataFileIF::DataType& a_dataType,
                          const IntVect&              a_num)
{
  // Maximum index IntVect
  IntVect a_max = a_num;
  a_max -= IntVect::Unit;

  // Box where the data is defined
  Box dataBox(IntVect::Zero,a_max);

  // Check data type to be read
  if (a_dataType == DataFileIF::ASCII)
  {
    // ASCII data - read the data one entry at a time

    // The data has one component
    RefCountedPtr<FArrayBox> data(new FArrayBox(dataBox,1));

    // Start all loop counters at zero
    IntVect state(IntVect::Zero);
    int comp = 0;

    bool first = true;
    a_maxValue = 0.0;

    // Loop until all loop counters are done
    while (state[SpaceDim-1] < a_num[SpaceDim-1])
    {
      Real curInput;

      if (a_dataType == DataFileIF::ASCII)
      {
        // Read the next ASCII value and store the data
        a_file >> curInput;
      }
      else
      if (a_dataType == DataFileIF::Binary)
      {
        // Read the next binary byte and store the data
        unsigned char curChar;
        a_file.read((char *)(&curChar),1);

        curInput = curChar;
      }

      (*data)(state,comp) = curInput;

      // Update the maximum data value
      if (first || a_maxValue < curInput)
      {
        first = false;
        a_maxValue = curInput;
      }

      // Increment the indices as needed - increasing x most rapidly
      for (int i = 0; i < SpaceDim; i++)
      {
        // Increment the counter
        state[i]++;

        // Is the counter in bounds?
        if (state[i] < a_num[i])
        {
          // If so, continue reading data
          break;
        }
        else
        {
          // If not, zero the counter and increment the next (nested loops)
          // unless this is the outermost loop counter
          if (i < SpaceDim-1)
          {
            state[i] = 0;
          }
        }
      }
    }

    m_ascii_data = data;
    m_binary_data = RefCountedPtr<BaseFab<unsigned char> >(NULL);
  }
  else if (a_dataType == DataFileIF::Binary)
  {
    // Binary data - read the data one entry at a time

    // The data has one component
    RefCountedPtr<BaseFab<unsigned char> > data(new BaseFab<unsigned char>(dataBox,1));

    // Start all loop counters at zero
    IntVect state(IntVect::Zero);
    int comp = 0;

    bool first = true;
    a_maxValue = 0.0;

    // Loop until all loop counters are done
    while (state[SpaceDim-1] < a_num[SpaceDim-1])
    {
      Real curInput;

      if (a_dataType == DataFileIF::ASCII)
      {
        // Read the next ASCII value and store the data
        a_file >> curInput;
      }
      else
      if (a_dataType == DataFileIF::Binary)
      {
        // Read the next binary byte and store the data
        unsigned char curChar;
        a_file.read((char *)(&curChar),1);

        curInput = curChar;
      }

      (*data)(state,comp) = curInput;

      // Update the maximum data value
      if (first || a_maxValue < curInput)
      {
        first = false;
        a_maxValue = curInput;
      }

      // Increment the indices as needed - increasing x most rapidly
      for (int i = 0; i < SpaceDim; i++)
      {
        // Increment the counter
        state[i]++;

        // Is the counter in bounds?
        if (state[i] < a_num[i])
        {
          // If so, continue reading data
          break;
        }
        else
        {
          // If not, zero the counter and increment the next (nested loops)
          // unless this is the outermost loop counter
          if (i < SpaceDim-1)
          {
            state[i] = 0;
          }
        }
      }
    }

    m_ascii_data = RefCountedPtr<FArrayBox>(NULL);
    m_binary_data = data;
  }
  else
  {
    // Unknown data - an error
    MayDay::Abort("Unknow data type specified for data file");
  }
}

void DataFileIF::MakeCorners(void)
{
  // Make an IntVectSet corresponding to the corners of a unit cube

  // Make the end points of a unit cube box
  IntVect lo(IntVect::Zero);
  IntVect hi(IntVect::Unit);

  // Make the unit cube box
  Box unitBox(lo,hi);

  // Save all its members (corners) as an IntVectSet
  m_corners.define(unitBox);
}

#include "NamespaceFooter.H"
