import biothings.hub.dataload.uploader as uploader
from hub.dataload.uploader import SnpeffPostUpdateUploader

class WellderlyFactoryUploader(uploader.DummySourceUploader,SnpeffPostUpdateUploader):
    """Data originally coming from: http://www.stsiweb.org/wellderly"""

    name = "wellderly"
    __metadata__ = {"mapper" : 'observed',
            "assembly" : "hg19",
            "src_meta" : {
                "url" : "https://genomics.scripps.edu/browser/",
                "license_url" : "https://genomics.scripps.edu/browser/page-help.html",
                "license_url_short": "https://goo.gl/e8OO17"
                }
            }

    #split_collections = ["wellderly_cg1","wellderly_cg10","wellderly_cg11",
    #                     "wellderly_cg12","wellderly_cg13","wellderly_cg14",
    #                     "wellderly_cg15","wellderly_cg16","wellderly_cg17",
    #                     "wellderly_cg18","wellderly_cg19","wellderly_cg2",
    #                     "wellderly_cg20","wellderly_cg21","wellderly_cg22",
    #                     "wellderly_cg3","wellderly_cg4","wellderly_cg5",
    #                     "wellderly_cg6","wellderly_cg7","wellderly_cg8",
    #                     "wellderly_cg9","wellderly_cgX","wellderly_cgY",
    #                     "wellderly_cgY1"]

    #@classmethod
    #def create(klass, db_conn_info, data_root, *args, **kwargs):
    #    return [klass(db_conn_info, data_root, collection_name=c,*args,*kwargs) for c in klass.split_collections]

    @classmethod
    def get_mapping(klass):
        mapping = {
            "wellderly": {
                "properties": {
                    "chrom": {
                        "type": "text",
                        "analyzer": "string_lowercase"
                    },
                    "pos": {
                        "type": "long"
                    },
                    "hg19": {
                        "properties": {
                            "start": {
                                "type": "integer"
                            },
                            "end": {
                                "type": "integer"
                            }
                        }
                    },
                    "ref": {
                        "type": "text",
                        "analyzer": "string_lowercase"
                    },
                    "alt": {
                        "type": "text",
                        "analyzer": "string_lowercase"
                    },
                    "vartype": {
                        "type": "text",
                        "analyzer": "string_lowercase"
                    },
                    # "alleles": {
                    #     "properties": {
                    #         "allele": {
                    #             "type": "text",
                    #             "analyzer": "string_lowercase"
                    #         },
                    #         "allele": {
                    #             "type": "float"
                    #         }
                    #     }
                    # },
                    "gene": {
                        "type": "text",
                        "analyzer": "string_lowercase",
                        "copy_to" : ["all"]
                    },
                    "coding_impact": {
                        "type": "text",
                        "analyzer": "string_lowercase"
                    },
                    "polyphen": {
                        "type": "text",
                        "analyzer": "string_lowercase"
                    },
                    "sift": {
                        "type": "text",
                        "analyzer": "string_lowercase"
                    }
                }
            }
        }
        return mapping

