import logging
import requests


class GraphqlContributorAvatars():

    OPERATION_NAME = "ContributorAvatars"

    QUERY = """
        query ContributorAvatars($subscribable: SubscribableInput!) {
        contributors(subscribable: $subscribable) {
            editors {
            ...ContributorFields
            __typename
            }
            curators {
            ...ContributorFields
            __typename
            }
            __typename
        }
        }

        fragment ContributorFields on ContributingUser {
        user {
            id
            profileImagePath(size: 12)
            __typename
        }
        uniqueActions {
            action
            count
            __typename
        }
        lastActionDate
        totalActionCount
        __typename
        }
    """

    def gql(self, variant_id: int):
        query = {
            "operationName": self.OPERATION_NAME,
            "query": self.QUERY,
            "variables": {
                "subscribable": {"id": variant_id, "entityType": "VARIANT"}
            }
        }
        return query

    def fetch(self, variant_id: int):
        try:
            response = requests.post(
                'https://civicdb.org/api/graphql',
                json=self.gql(variant_id=variant_id)
            )
            return response.json()
        except Exception as e:
            logging.error(f"Error in {self.OPERATION_NAME}: {e}")
            raise
