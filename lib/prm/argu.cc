#include "argu.hh"

using namespace prm;

tag::Base * tag::Key::dup() const
{
	return new Key(key);
}

tag::Key::Key(int a_key)
{
	key = a_key;
}

argu::Base::Base(const std::string & n, TagSable & t) :
	name(n),
	ts(t),
	set_altering(0)
{
}

void argu::Base::dump(std::ostream & output) const
{
	if (const tag::Name * n = ts.find<tag::Name>()) output << n->text;
	else output << name;
	output << '=' << ts.to_str();
}

void argu::Base::set(const std::string & str)
{
	ts.read_str(str);
}

void argu::Base::addto_parser(arg::Parser & parser)
{
	int key = 0;
	if (const tag::Key * k = ts.find<tag::Key>()) key = k->key;
	arg::Option & opt = parser.add_opt(key, get_name())
		.set(& set_altering, true)
		.show_default();
	if (const tag::Desc * d = ts.find<tag::Desc>()) opt.help(d->text);
	addto_option(opt);
}

bool argu::Base::alter()
{
	if (set_altering) do_altering();
	return set_altering;
}

template <> std::string argu::type_word<bool>() {return "BOOL";}
template <> std::string argu::type_word<int>() {return "INT";}
template <> std::string argu::type_word<size_t>() {return "SIZE";}
template <> std::string argu::type_word<double>() {return "DOUBLE";}

Argu::~Argu()
{
	clear();
}

void Argu::add(Param & p)
{
	if (argus.size()) argus.push_back(0);
	for (std::vector<Param::VarEntry>::iterator i = p.vars.begin(); i != p.vars.end(); i ++) if (i->var->find<tag::CmdLine>()) {
		argu::Base * c = 0;
		// The Type infomation for Var<Type> can not be recovered.
		// So, we are making tabulated translations:
		if (Var<bool> * v = dynamic_cast<Var<bool> *>(i->var)) c = new argu::Arg<bool>(i->name, * v);
		else if (Var<int> * v = dynamic_cast<Var<int> *>(i->var)) c = new argu::Arg<int>(i->name, * v);
		else if (Var<size_t> * v = dynamic_cast<Var<size_t> *>(i->var)) c = new argu::Arg<size_t>(i->name, * v);
		else if (Var<double> * v = dynamic_cast<Var<double> *>(i->var)) c = new argu::Arg<double>(i->name, * v);
		else if (Var<std::string> * v = dynamic_cast<Var<std::string> *>(i->var)) c = new argu::Arg<std::string>(i->name, * v);
		if (c) argus.push_back(c);
	}
}

void Argu::clear()
{
	// This invalidates the SubParser created by make_parser
	while (argus.size()) {
		if (argu::Base * b = argus.back()) delete b;
		argus.pop_back();
	}
}

arg::SubParser * Argu::make_parser(std::string title)
{
	arg::SubParser * sp = new arg::SubParser;
	if (title.size()) sp->add_help(title + '\n');
	sp->add_help("Available parameters:");
	sp->add_help("");
	for (std::vector<argu::Base *>::iterator i = argus.begin(); i != argus.end(); i ++) {
		if (argu::Base * b = * i) b->addto_parser(* sp);
		else sp->add_help("");
	}
	sp->add_help("");
	sp->add_opt_help();
	return sp;
}

void Argu::addto_parser(arg::Parser & parser)
{
	for (std::vector<argu::Base *>::iterator i = argus.begin(); i != argus.end(); i ++) {
		if (argu::Base * b = * i) b->addto_parser(parser);
		else parser.add_help("");
	}
}

bool Argu::alter()
{
	bool any_change = false;
	for (std::vector<argu::Base *>::iterator i = argus.begin(); i != argus.end(); i ++) {
		if (argu::Base * b = * i) any_change = b->alter() || any_change;
	}
	return any_change;
}

size_t Argu::size() const
{
	return argus.size();
}
