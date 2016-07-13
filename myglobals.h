#pragma once
//#include <vld.h>
#include <memory>
#include <sstream>
#include <iostream>

#include <ctime>
#include <QTCore/QtConcurrentRun>
#include <QTCore/QString>
#include <QFuture>

#include "Libraries/Logger.h"

#define qtawait(theType, theCall) ([&]()->theType{	\
		QEventLoop el;	\
		QFuture<theType> theResult;\
		auto threadLambda = [&]()->theType{			\
			bool succ = QObject::connect(QThread::currentThread(), SIGNAL(finished()), &el, SLOT(quit()));	\
			theType result = theCall;							\
			el.quit();								\
			qDebug() << "finished async function with return value" ;	\
			return result;										\
		};		\
		theResult = QtConcurrent::run((function<theType()>)threadLambda);		\
		el.exec();									\
		return theResult.result();					\
	}())

#define qtawait_void(theCall) ([&]{	\
		QEventLoop el;						\
		auto threadLambda = [&]{			\
			bool succ = QObject::connect(QThread::currentThread(), SIGNAL(finished()), &el, SLOT(quit()));	\
			theCall;								\
			el.quit();								\
			qDebug() << "finished void async function";\
		};		\
		QtConcurrent::run((function<void()>)threadLambda);		\
		el.exec();									\
	}())

#define qtawait_local_statments(theBody) {	\
		QEventLoop el;						\
		auto threadLambda = [&]{			\
			try								\
			{								\
				bool succ = QObject::connect(QThread::currentThread(), SIGNAL(finished()), &el, SLOT(quit()));	\
				theBody;								\
				el.quit();								\
				/*qDebug() << "finished async local statments";*/	\
			}	\
			catch(std::exception e) \
			{\
				std::cout<<"Exception in qtawait_local_statements : "<<e.what() <<endl; \
				throw e;											\
			}		\
		};		\
		QtConcurrent::run(threadLambda);		\
		/*QtConcurrent::run((function<void()>)threadLambda);		*/\
		el.exec();									\
	}


#ifdef LOG_ENABLED
	#define stream_to_log(theStream) {stream_to_log_inner(theStream)}
	#define stream_to_log_sameLine(theStream) {stream_to_log_inner_sameLine(theStream)}
#else
	#define stream_to_log(theStream)
	#define stream_to_log_sameLine(theStream)
#endif

#define stream_to_log_inner(theStream) { using namespace HyperSpectralToolbox; \
					auto sp_ofstream = Logger::openStream(); \
					(*sp_ofstream) << theStream << "\n";	\
					sp_ofstream->close(); \
					std::cout << theStream <<endl; \
					}					
			
#define stream_to_log_inner_sameLine(theStream) { using namespace HyperSpectralToolbox; \
					auto sp_ofstream = Logger::openStream(); \
					(*sp_ofstream) << theStream;	\
					sp_ofstream->close(); \
					std::cout << theStream; \
					}	

//#define ts(theTokens) ([&]()->std::string { std::stringstream sst; sst << theTokens; return sst.str();}())

template <typename T> 
inline std::string ts( const T & theToken ) { 
	std::ostringstream os; 
	os << theToken; 
	return os.str(); 
} 

template <typename T> 
inline QString tqs( const T & theToken ) { 
	std::ostringstream os; 
	os << theToken; 
	return QString(os.str().c_str()); 
} 

#define shnew(theType, constructorArgs) \
	make_shared<theType>(theType constructorArgs )

#define shptr(theType) std::shared_ptr<theType>

#include "Libraries/HSWO.h"

struct MyFutureResult{
	std::shared_ptr<HyperSpectralToolbox::HSWO> spHswo;
	QFuture<void> future;
};