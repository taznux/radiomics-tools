#include "gdcmReader.h"
#include "gdcmAttribute.h"

#include "gdcmGlobal.h"
#include "gdcmDicts.h"
#include "gdcmDict.h"
#include "gdcmCSAHeader.h"
#include "gdcmPrivateTag.h"



void recursive_read(const gdcm::DataSet& ds)
{
    static int depth = 0;
    const gdcm::Global& g = gdcm::Global::GetInstance(); // sum of all knowledge !
    const gdcm::Dicts &dicts = g.GetDicts();
    const gdcm::Dict &dict = dicts.GetPublicDict();

    for (gdcm::DataSet::ConstIterator it = ds.Begin(); it!=ds.End(); ++it)
    {
        const gdcm::DataElement& elem = *it;

        if(!elem.IsEmpty())
        {    
            const gdcm::Tag& tag = elem.GetTag();

            for(int d = 0; d < depth; d++)
                std::cout << ' ';
            std::cout << dict.GetDictEntry(tag).GetKeyword() << ": ";
            std::cout << elem.GetValue() << std::endl;
			try
			{
				gdcm::SmartPointer<gdcm::SequenceOfItems> ssqi = elem.GetValueAsSQ();
				if (ssqi != NULL)
				{
					for (int i = 1; i <= ssqi->GetNumberOfItems(); i++)
					{
						const gdcm::Item & item = ssqi->GetItem(i); // Item start at #1
						const gdcm::DataSet& nested_ds = item.GetNestedDataSet();

						depth++;
						recursive_read(nested_ds);
						depth--;
					}
				}
            }
            catch (const gdcm::Exception& e)
            {
                std::cerr << e.what() << std::endl;
            }
        }
     }
}

int main( int argc, char* argv[] )
{
    if( argc < 2 )
    {
        std::cerr << "Usage: " << argv[0] << " DicomFilename " << std::endl;
        return EXIT_FAILURE;
    }

    gdcm::Reader reader;
    reader.SetFileName(argv[1]);
    if( !reader.Read() )
    {
        std::cerr << "Cannot read DICOM file!!" << std::endl;
        return 1;
    }

    gdcm::MediaStorage ms;
    ms.SetFromFile( reader.GetFile() );
    if( ms != gdcm::MediaStorage::SpacialRegistrationStorage)
    {
        std::cerr << "No DICOM REG!!" << std::endl;
    }

    const gdcm::DataSet& ds = reader.GetFile().GetDataSet();
    recursive_read(ds);
    

    return EXIT_SUCCESS;
}
