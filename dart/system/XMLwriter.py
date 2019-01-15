#!/usr/bin/python2.7

"""
==========================================================================================

Author:				By Bruce Eckel, Permission is granted to use or modify without 
					payment as long as this copyright notice is retained.
Copyright (C):		2006 MindView Inc. www.MindView.net
DART version:		1.2 (25-11-2008)
DART module: 		XMLwriter.py
Module function:	Rapidly assemble XML using minimal coding.

					Everything is a Node, and each Node can either have a value 
					or subnodes. Subnodes can be appended to Nodes using '+=', 
					and a group of Nodes can be strung together using '+'.

					Create a node containing a value by saying 
					Node("tag", "value")
					You can also give attributes to the node in the constructor:
					Node("tag", "value", attr1 = "attr1", attr2 = "attr2")
					or without a value:
					Node("tag", attr1 = "attr1", attr2 = "attr2")

					To produce xml from a finished Node n, say n.xml() (for 
					nicely formatted output) or n.rawxml().

					You can read and modify the attributes of an xml Node using 
					getAttribute(), setAttribute(), or delAttribute().

					You can find the value of the first subnode with tag == "tag"
					by saying n["tag"]. If there are multiple instances of n["tag"],
					this will only find the first one, so you should use node() or
					nodeList() to narrow your search down to a Node that only has
					one instance of n["tag"] first.

					You can replace the value of the first subnode with tag == "tag"
					by saying n["tag"] = newValue. The same issues exist as noted
					in the above paragraph.

					You can find the first node with tag == "tag" by saying 
					node("tag"). If there are multiple nodes with the same tag 
					at the same level, use nodeList("tag").

					The Node class is also designed to create a kind of "domain 
					specific language" by subclassing Node to create Node types 
					specific to your problem domain.

					This implementation uses xml.dom.minidom which is available
					in the standard Python 2.4 library. However, it can be 
					retargeted to use other XML libraries without much effort.
Module depenencies:	python2.4 xml.dom.minidom modules

for further information, please contact:
			- Bruce Eckel, MindView Inc. www.MindView.net

==========================================================================================
"""

"""Import modules"""
from xml.dom.minidom import getDOMImplementation, parseString
import copy, re

class Node(object):
    """
    Everything is a Node. The XML is maintained as (very efficient)
    Python objects until an XML representation is needed.
    """
    def __init__(self, tag, value = None, **attributes):
        self.tag = tag.strip()
        self.attributes = attributes
        self.children = []
        self.value = value
        if self.value:
            self.value = self.value.strip()

    def getAttribute(self, name):
        """
        Read XML attribute of this node.
        """
        return self.attributes[name]

    def setAttribute(self, name, item):
        """
        Modify XML attribute of this node.
        """
        self.attributes[name] = item

    def delAttribute(self, name):
        """
        Remove XML attribute with this name.
        """
        del self.attributes[name]

    def node(self, tag):
        """ 
        Recursively find the first subnode with this tag. 
        """
        if self.tag == tag:
            return self
        for child in self.children:
            result = child.node(tag)
            if result:
                return result
        return False
        
    def nodeList(self, tag):
        """ 
        Produce a list of subnodes with the same tag. 
        Note:
        It only makes sense to do this for the immediate 
        children of a node. If you went another level down, 
        the results would be ambiguous, so the user must 
        choose the node to iterate over.
        """
        return [n for n in self.children if n.tag == tag]

    def __getitem__(self, tag):
        """ 
        Produce the value of a single subnode using operator[].
        Recursively find the first subnode with this tag. 
        If you want duplicate subnodes with this tag, use
        nodeList().
        """
        subnode = self.node(tag)
        if not subnode:
            raise KeyError
        return subnode.value

    def __setitem__(self, tag, newValue):
        """ 
        Replace the value of the first subnode containing "tag"
        with a new value, using operator[].
        """
        assert isinstance(newValue, str), "Value " + str(newValue) + " must be a string"
        subnode = self.node(tag)
        if not subnode:
            raise KeyError
        subnode.value = newValue

    def __iadd__(self, other):
        """
        Add child nodes using operator +=
        """
        assert isinstance(other, Node), "Tried to += " + str(other)
        self.children.append(other)
        return self

    def __add__(self, other):
        """
        Allow operator + to combine children
        """
        return self.__iadd__(other)

    def __str__(self):
        """
        Display this object (for debugging)
        """
        result = self.tag + "\n"
        for k, v in self.attributes.items():
            result += "    attribute: %s = %s\n" % (k, v)
        if self.value:
            result += "    value: [%s]" % self.value
        return result
        
    # The following are the only methods that rely on the underlying
    # Implementation, and thus the only methods that need to change
    # in order to retarget to a different underlying implementation.

    # A static dom implementation object, used to create elements:        
    doc = getDOMImplementation().createDocument(None, None, None)

    def dom(self):
        """
        Lazily create a minidom from the information stored
        in this Node object.
        """
        element = Node.doc.createElement(self.tag)
        for key, val in self.attributes.items():
            element.setAttribute(key, val)
        if self.value:
            assert not self.children, "cannot have value and children: " + str(self)
            element.appendChild(Node.doc.createTextNode(self.value))
        else:
            for child in self.children:
                element.appendChild(child.dom()) # Generate children as well
        return element

    def xml(self, separator = '  '):
        return self.dom().toprettyxml(separator)

    def rawxml(self):
        return self.dom().toxml()

    #staticmethod
    def create(dom):
        """
        Create a Node representation, given either
        a string representation of an XML doc, or a dom.
        """
        if isinstance(dom, str):
            # Strip all extraneous whitespace so that
            # text input is handled consistently:
            dom = re.sub("\s+", " ", dom)
            dom = dom.replace("> ", ">")
            dom = dom.replace(" <", "<")
            return Node.create(parseString(dom))
        if dom.nodeType == dom.DOCUMENT_NODE:
            return Node.create(dom.childNodes[0])
        if dom.nodeName == "#text":
            return
        node = Node(dom.nodeName)
        if dom.attributes:
            for name, val in dom.attributes.items():
                node.setAttribute(name, val)
        for n in dom.childNodes:
            if n.nodeType == n.TEXT_NODE and n.wholeText.strip():
                node.value = n.wholeText
            else:
                subnode = Node.create(n)
                if subnode:
                    node += subnode
        return node
