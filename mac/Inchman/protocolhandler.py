#!/usr/bin/python
#
# Inchman Protocol Handler
#
# Created by Aidan Lane, 23 Jan 2013
# Copyright (c) 2013 Monash University
#
# Note that this script download the xml file into the current working directory.

import sys, re, urllib2
import paramsweeper


def id_is_valid(id):
    uuid4hex = re.compile('[0-9a-f]{32}\Z', re.I)
    return uuid4hex.match(id) != None


def fetch_xml_file(inchman_exec, host, id):
    url = 'http://' + host + '/api/modelSaveDownload?id=' + id
    print "Trying %s ..." % host
    response = urllib2.urlopen(url)
    filename = id + '.xml'
    f = open(filename, 'w')
    f.write(response.read())
    f.close()
    response.close()
    return filename


def simulate(inchman_exec, id):
    # Try local server first
    # IMPORTANT use hard-coded host so that the 'inchman://' handler isn't exploited to load files from other locations
    filename = None
    try:
        filename = fetch_xml_file(inchman_exec, 'localhost:8888', id)
    except urllib2.URLError, urllib2.HTTPError:
        try:
            # Now try remote server
            # IMPORTANT use hard-coded host so that the 'inchman://' handler isn't exploited to load files from other locations
            filename = fetch_xml_file(inchman_exec, 'inchman-editor.appspot.com', id)
        except urllib2.URLError, urllib2.HTTPError:
            raise Exception("Failed to find the experiment on any server")

    paramsweeper.exec_sweep(inchman_exec, filename)


def handle_url(inchman_exec, complete_url):
    handler, fullPath = complete_url.split(':', 1)

    # Validate handler
    if handler != 'inchman':
        raise Exception("Invalid handler: %s" % handler)

    action, id = fullPath.strip('/').split('/')

    # Validate action
    if action != 'simulate':
        raise Exception("Unknown action %s" % action)

    # Validate id - ensure that the id doesn't contain any malicious characters
    if not id_is_valid(id):
        raise Exception("Invalid experiment ID: %s" % id)

    simulate(inchman_exec, id)


def main(argv):
    handle_url(argv[1], argv[2])


if __name__ == '__main__':
    main(sys.argv)
