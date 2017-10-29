{
  "baseCommand": [
    "/opt/bin/python",
    "/opt/MatchedFilter/scanMatchedFilter.py"
  ],
  "sbg:cmdPreview": "/opt/bin/python /opt/MatchedFilter/scanMatchedFilter.py  /opt/MatchedFilter/Training/S2_H3K27ac_positivesAll_MFscores.bed  /opt/MatchedFilter/Training/S2_H3K27ac_negativesAll_MFscores.bed",
  "id": "https://api.sbgenomics.com/v2/apps/anurag.sethi/matchedfilter/scanmatchedfilter/3/raw/",
  "sbg:projectName": "MatchedFilter",
  "sbg:image_url": null,
  "sbg:validationErrors": [],
  "sbg:latestRevision": 3,
  "successCodes": [],
  "outputs": [
    {
      "type": [
        "null",
        {
          "type": "array",
          "items": "File"
        }
      ],
      "id": "#SVM_predictions",
      "label": "SVM_preditions",
      "outputBinding": {
        "glob": "*SVM*",
        "sbg:inheritMetadataFrom": "#fileList"
      }
    }
  ],
  "label": "scanMatchedFilter",
  "sbg:id": "anurag.sethi/matchedfilter/scanmatchedfilter/3",
  "sbg:job": {
    "inputs": {
      "outputPrefix": "outputPrefix-string-value",
      "metaProfileList": {
        "size": 0,
        "secondaryFiles": [],
        "path": "/path/to/metaProfileList.ext",
        "class": "File"
      },
      "peakFileList": {
        "size": 0,
        "secondaryFiles": [],
        "path": "/path/to/peakFileList.ext",
        "class": "File"
      },
      "chrNameFile": {
        "size": 0,
        "secondaryFiles": [],
        "path": "/path/to/chrNameFile.ext",
        "class": "File"
      },
      "fileList": {
        "size": 0,
        "secondaryFiles": [],
        "path": "/path/to/fileList.ext",
        "class": "File"
      }
    },
    "allocatedResources": {
      "cpu": 1,
      "mem": 8000
    }
  },
  "sbg:modifiedBy": "anurag.sethi",
  "stdout": "",
  "stdin": "",
  "requirements": [],
  "temporaryFailCodes": [],
  "sbg:revisionsInfo": [
    {
      "sbg:revision": 0,
      "sbg:modifiedOn": 1509276304,
      "sbg:revisionNotes": null,
      "sbg:modifiedBy": "anurag.sethi"
    },
    {
      "sbg:revision": 1,
      "sbg:modifiedOn": 1509276803,
      "sbg:revisionNotes": "Added inputs",
      "sbg:modifiedBy": "anurag.sethi"
    },
    {
      "sbg:revision": 2,
      "sbg:modifiedOn": 1509279893,
      "sbg:revisionNotes": null,
      "sbg:modifiedBy": "anurag.sethi"
    },
    {
      "sbg:revision": 3,
      "sbg:modifiedOn": 1509282746,
      "sbg:revisionNotes": null,
      "sbg:modifiedBy": "anurag.sethi"
    }
  ],
  "hints": [
    {
      "value": 1,
      "class": "sbg:CPURequirement"
    },
    {
      "value": 8000,
      "class": "sbg:MemRequirement"
    },
    {
      "dockerPull": "images.sbgenomics.com/anurag_sethi/matchedfilter:1",
      "dockerImageId": "",
      "class": "DockerRequirement"
    }
  ],
  "inputs": [
    {
      "type": [
        "null",
        "File"
      ],
      "inputBinding": {
        "separate": true,
        "sbg:cmdInclude": true,
        "position": 1
      },
      "label": "file list",
      "id": "#fileList",
      "sbg:stageInput": null,
      "description": "file with the list of chromatin signals in the format (2 column tab delimited)"
    },
    {
      "type": [
        "null",
        "File"
      ],
      "inputBinding": {
        "separate": true,
        "sbg:cmdInclude": true,
        "position": 2
      },
      "label": "metaprofile list",
      "id": "#metaProfileList",
      "sbg:stageInput": null,
      "description": "File with list of metaprofiles in the format (2 column tab delimited) containing signaltype and filename"
    },
    {
      "type": [
        "null",
        "File"
      ],
      "inputBinding": {
        "separate": true,
        "sbg:cmdInclude": true,
        "position": 3
      },
      "label": "chromosome name file",
      "id": "#chrNameFile",
      "sbg:stageInput": null,
      "description": "chromosome size file"
    },
    {
      "type": [
        "null",
        "File"
      ],
      "inputBinding": {
        "separate": true,
        "sbg:cmdInclude": true,
        "position": 4
      },
      "label": "peak file list",
      "id": "#peakFileList",
      "sbg:stageInput": null,
      "description": "name of the peak file list for each signal type"
    },
    {
      "type": [
        "null",
        "string"
      ],
      "inputBinding": {
        "separate": true,
        "sbg:cmdInclude": true,
        "position": 7
      },
      "label": "output prefix",
      "id": "#outputPrefix",
      "sbg:stageInput": null,
      "description": "prefix for all output files"
    }
  ],
  "sbg:publisher": "sbg",
  "sbg:revision": 3,
  "sbg:modifiedOn": 1509282746,
  "arguments": [
    {
      "separate": true,
      "valueFrom": "/opt/MatchedFilter/Training/S2_H3K27ac_positivesAll_MFscores.bed",
      "prefix": "",
      "position": 5
    },
    {
      "separate": true,
      "valueFrom": "/opt/MatchedFilter/Training/S2_H3K27ac_negativesAll_MFscores.bed",
      "position": 6
    }
  ],
  "sbg:contributors": [
    "anurag.sethi"
  ],
  "sbg:createdBy": "anurag.sethi",
  "sbg:project": "anurag.sethi/matchedfilter",
  "class": "CommandLineTool",
  "cwlVersion": "sbg:draft-2",
  "sbg:appVersion": [
    "sbg:draft-2"
  ],
  "sbg:sbgMaintained": false,
  "description": "",
  "sbg:createdOn": 1509276304
}
