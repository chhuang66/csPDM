#ifndef CTRL_TAG_HH
#define CTRL_TAG_HH 1
#include "prm.hh"
#include "name_set.hh"
namespace prm {
	namespace tag {
		class Control : // value can be interactively controlled
			virtual public Base
		{
			Base * dup() const {return new Control;}
		};

		class Step : // a numerical value
			virtual public Control
		{
			Base * dup() const {return new Step(step, page);}
		public:
			const double step;
			const double page;
			Step(double st, double pg) : step(st), page(pg) {}
		};

		class Select :
			virtual public Control
		{
			Base * dup() const {return new Select(names);}
		public:
			NameSet & names;
			Select(NameSet & ns) : names(ns) {}
		};
	}
}
#endif // CTRL_TAG_HH
