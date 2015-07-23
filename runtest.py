#!/usr/bin/env python
import unittest
import optparse
import sys

def get_option_parser():
    parser = optparse.OptionParser()
    parser.add_option('-v', '--verbosity',
        dest='verbosity',
        default=1,
        type='int',
        help='Verbosity of output')

    parser.add_option("-a","--all",
        action="store_true", 
        dest="test_all",
        default=False,
        help="Run all tests")

    parser.add_option("-m","--math",
        action="store_true", 
        dest="test_math",
        default=False,
        help="Run tests from math engine")

    parser.add_option("-f","--flask",
        action="store_true", 
        dest="test_flask",
        default=False,
        help="Run tests for flask")
    return parser

def load_modules(options):
    modules = []
    if options.test_math:
        from callculation.tests import test_frame
        modules.append(test_frame)
    if options.test_flask:
        from tests import test_app
        modules.append(test_app)
    return modules


def runtests(suite, verbosity=1):
    results = unittest.TextTestRunner(verbosity=verbosity).run(suite)
    return results.failures, results.errors

if __name__ == '__main__':
    parser = get_option_parser()
    options, args = parser.parse_args()

    suite = unittest.TestSuite()
    for module in load_modules(options):
        print 'Adding tests for "%s"' % module.__name__
        module_suite = unittest.TestLoader().loadTestsFromModule(module)
        suite.addTest(module_suite)

    failures, errors = runtests(suite, options.verbosity)

    if errors:
        sys.exit(2)
    elif failures:
        sys.exit(1)

    sys.exit(0)