#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <memory>

#include "definitions.h"

#ifdef LOG_ENABLED
	#define LOG_START if(true) {
#else
	#define LOG_START if(false) {
#endif

#define LOG_END }

#define LOG_DISABLE if(false) {

namespace HyperSpectralToolbox
{	
	class Logger
	{
	protected:
		//static std::ofstream out_stream;		

	public:
		static void InitializeLogger()
		{	
			std::ofstream out_stream;		
			out_stream.open(getLogFileName().c_str());			
			out_stream.close();
		}
		
		static std::string getLogFileName()
		{
			using namespace std;
			return string(logFileFolder) + string("\\") +string(logFileName);
		}

		static std::shared_ptr<std::ofstream> openStream()
		{
			std::shared_ptr<std::ofstream> sp_out_stream(new std::ofstream(getLogFileName().c_str(), std::ios::app));
			//if(out_stream.is_open() != true) out_stream.open(getLogFileName().c_str());			
			return sp_out_stream;
		}		

		static void closeStream()
		{
			//out_stream.close();
		}

		inline static void stringLine(char* line)
		{
			LOG_START
				std::shared_ptr<std::ofstream> sp_out_stream(new std::ofstream(getLogFileName().c_str()));
				(*sp_out_stream) << line<<"\n";			
			LOG_END
		}

		inline static void sstring(char* str)
		{
			LOG_START
				std::shared_ptr<std::ofstream> sp_out_stream(new std::ofstream(getLogFileName().c_str()));
				(*sp_out_stream)<<str;			
			LOG_END
		}

		inline static void decimal(double line)
		{
			LOG_START
				std::shared_ptr<std::ofstream> sp_out_stream(new std::ofstream(getLogFileName().c_str()));
				(*sp_out_stream)<<line;			
			LOG_END
		}

		inline static void decimal(double line, int width)
		{
			LOG_START
				std::shared_ptr<std::ofstream> sp_out_stream(new std::ofstream(getLogFileName().c_str()));
				(*sp_out_stream).width(width);
				(*sp_out_stream)<<line;				
			LOG_END
		}

		inline static void endl()
		{
			LOG_START
				std::shared_ptr<std::ofstream> sp_out_stream(new std::ofstream(getLogFileName().c_str()));
				(*sp_out_stream)<<"\n";			
			LOG_END
		}		
	};	
}