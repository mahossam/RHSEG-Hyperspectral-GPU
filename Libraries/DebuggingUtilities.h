#pragma once

#include <string>
#include <map>
#include <vector>
#include <set>

namespace HyperSpectralToolbox
{

	template <typename ObjectType>
	class DebuggingUtilities
	{		
		//static map<string, int> variablesInt_;
		//static map<string, double> variablesDecimal_;
		typedef ObjectType DebuggingUtilities_ObjectType;
		typedef std::set<ObjectType>* TSetPtr;
		typedef std::map<std::string, TSetPtr> mapStringToSet;
		
		std::map<std::string, std::set<DebuggingUtilities_ObjectType>*> variablesArray_;
		DebuggingUtilities(){}

	public:		
		static DebuggingUtilities<DebuggingUtilities_ObjectType>& getInstance()
		{
			static DebuggingUtilities<DebuggingUtilities_ObjectType> singleInstance;
			return singleInstance;
		}
		/*
		static void AddInteger(string name, int int_value)
		{
			variablesInt_.insert(make_pair(name, int_value));
		}
		static void AddDecimal(string name, double decimal_value)
		{
			variablesDecimal_.insert(make_pair(name, decimal_value));
		}
		*/

		void AddArray(std::string name)
		{
			variablesArray_.insert(make_pair(name, new std::set<DebuggingUtilities_ObjectType>()));
		}

		void AddArrayValue(std::string name, DebuggingUtilities_ObjectType object_value)
		{
			std::set<DebuggingUtilities_ObjectType>* pArray=NULL;
			/*
			if( (pArray = GetArray(name)) != NULL)
			{
				(*pArray).insert(object_value);
			}
			*/
		}

/*
		std::set<DebuggingUtilities_ObjectType>* GetArray(std::string name)
		{
			//std::map<string, set<DebuggingUtilities_ObjectType>*>::iterator it = variablesArray_.find(name);						
			mapStringToSet::iterator it = variablesArray_.find(name);
			if(it != variablesArray_.end()) return (it->second);
			return NULL;
	}
*/

		virtual ~DebuggingUtilities()
		{
			/*
			std::map<std::string, set<DebuggingUtilities_ObjectType>*>::reverse_iterator rev_it ;
			for ( rev_it = variablesArray_.rbegin(); rev_it  != variablesArray_.rend(); rev_it ++)
			{
				delete (DebuggingUtilities_ObjectType*)(rev_it->second);
				variablesArray_.erase(rev_it->first);
			}
			*/
		}


		/*
		static int GetInteger(string name)
		{
			map<string, int>::iterator it = variablesInt_.find(name);
			if(it != variablesInt_.end()) return *it;
			return INT_MAX;
		}
		static double GetDecimal(string name)
		{
			map<string, double>::iterator it = variablesDecimal_.find(name);
			if(it != variablesDecimal_.end()) return *it;
			return INT_MAX;
		}
		*/
	};
	
	typedef DebuggingUtilities<int> IntegerDebuggingUtilities;
}

